% NIfTIMRSRead
%   Code borrowed from Osprey's io_loadspec_niimrs.m, created by Georg
%   Oeltzschner, Johns Hopkins University (2021)
%
%   Adapted for Gannet by Mark Mikkelsen, Weill Cornell Medicine (2022)
%
%   Description:
%   Reads in MRS data stored according to the NIfTI-MRS format
%   (https://doi.org/10.1002/mrm.29418)
%
%   Dependencies:
%   This function requires the dcm2nii toolbox (Xiangrui Li) in the main
%   Gannet folder to be on the MATLAB search path
%   (https://github.com/xiangruili/dicm2nii)

function MRS_struct = NIfTIMRSRead(MRS_struct, fname, fname_w)

ii = MRS_struct.ii;

% Read in the data using the dicm2nii toolbox
% (https://github.com/xiangruili/dicm2nii)
try
    nii = nii_tool('load', fname);
catch ME
    switch ME.identifier
        case 'MATLAB:UndefinedFunction'
            error(['Cannot find the function ''nii_tool.m''. ' ...
                   'Please ensure that you have added the dcm2nii ', ...
                   'folder in the main Gannet folder to your MATLAB ', ...
                   'search path.']);
        otherwise
            rethrow(ME);
    end
end

% Extract the header and header extensions
hdr     = nii.hdr;
hdr_ext = jsondecode(nii.ext.edata_decoded);

% Extract the raw time-domain data
fids = double(nii.img);

MRS_struct.p.LarmorFreq(ii) = hdr_ext.SpectrometerFrequency;
MRS_struct.p.sw(ii)         = 1/hdr.pixdim(5);
if strcmp(hdr_ext.Manufacturer,'GE')
    MRS_struct.p.TE(ii)     = hdr_ext.EchoTime / 1e3;
    MRS_struct.p.TR(ii)     = hdr_ext.RepetitionTime / 1e3;
else
    MRS_struct.p.TE(ii)     = hdr_ext.EchoTime * 1e3;
    MRS_struct.p.TR(ii)     = hdr_ext.RepetitionTime * 1e3;
end
MRS_struct.p.voxdim(ii,:)   = hdr.pixdim(2:4);

% Specify dimensions
[dims, fids, npoints] = specify_dims(hdr, hdr_ext, fids);

sz = size(fids);

MRS_struct.p.npoints(ii) = npoints;
if dims.averages && dims.subSpecs
    MRS_struct.p.Navg(ii) = sz(dims.averages) * sz(dims.subSpecs);
elseif length(sz) == 2 && sz(2) == 1
    MRS_struct.p.Navg(ii) = 1;
else
    MRS_struct.p.Navg(ii) = sz(dims.averages);
end

if dims.averages && dims.subSpecs
    fids = reshape(fids, [sz(dims.t) sz(dims.coils) sz(dims.averages) * sz(dims.subSpecs)]);
    fids = permute(fids, [2 1 3]);
end

% Load water reference
if nargin == 3

    nii_w     = nii_tool('load', fname_w);
    hdr_w     = nii.hdr;
    hdr_w_ext = jsondecode(nii.ext.edata_decoded);

    fids_w = double(nii_w.img);

    MRS_struct.p.sw_water(ii)     = 1/hdr_w.pixdim(5);
    if strcmp(hdr_w_ext.Manufacturer,'GE')
        MRS_struct.p.TE_water(ii) = hdr_w_ext.EchoTime / 1e3;
        MRS_struct.p.TR_water(ii) = hdr_w_ext.RepetitionTime / 1e3;
    else
        MRS_struct.p.TE_water(ii) = hdr_w_ext.EchoTime * 1e3;
        MRS_struct.p.TR_water(ii) = hdr_w_ext.RepetitionTime * 1e3;
    end

    % Specify dimensions
    [dims_w, fids_w] = specify_dims(hdr_w, hdr_w_ext, fids_w);

    sz_w = size(fids_w);

    MRS_struct.p.npoints_water(ii) = sz_w(dims_w.t);
    if dims_w.averages && dims_w.subSpecs
        MRS_struct.p.Nwateravg(ii) = sz_w(dims_w.averages) * sz_w(dims_w.subSpecs);
    elseif length(sz_w) == 2 && sz_w(2) == 1
        MRS_struct.p.Nwateravg(ii) = 1;
    else
        MRS_struct.p.Nwateravg(ii) = sz_w(dims_w.averages);
    end

    if dims_w.averages && dims_w.subSpecs
        fids_w = reshape(fids_w, [sz_w(dims_w.t) sz_w(dims_w.coils) sz_w(dims_w.averages) * sz_w(dims_w.subSpecs)]);
        fids_w = permute(fids_w, [2 1 3]);
    end

end

if dims.coils > 0

    if nargin == 3
        [nCh, nPts, nReps] = size(fids_w);
        noise_pts = false(1,nPts);
        noise_pts(ceil(0.75*nPts):end) = true;
        noise_pts = repmat(noise_pts, [1 nReps]);
        tmp_fids_w = reshape(fids_w, [nCh nPts*nReps]);

        e = tmp_fids_w(:,noise_pts);
        Psi = e*e';
        fids_w_avg = mean(fids_w,3);
        S = fids_w_avg(:,1);
        w = (S'*(Psi\S))^-1 * S' / Psi;
        fids_w = w.' .* fids_w;
        MRS_struct.fids.data_water = mean(squeeze(sum(fids_w,1)),2);
    end

    [nCh, nPts, nReps] = size(fids);
    noise_pts = false(1,nPts);
    noise_pts(ceil(0.75*nPts):end) = true;
    noise_pts = repmat(noise_pts, [1 nReps]);
    tmp_fids = reshape(fids, [nCh nPts*nReps]);

    e = tmp_fids(:,noise_pts);
    Psi = e*e';
    if nargin == 2
        fids_avg = mean(fids,3);
        S = fids_avg(:,1);
    end
    w = (S'*(Psi\S))^-1 * S' / Psi;
    fids = w.' .* fids;
    MRS_struct.fids.data = squeeze(sum(fids,1));

else

    MRS_struct.fids.data = fids;
    if nargin == 3
        if dims_w.averages > 0
            MRS_struct.fids.data_water = mean(fids_w,2);
        else
            MRS_struct.fids.data_water = fids_w;
        end
    end

end

switch hdr_ext.Manufacturer

    case 'GE'
        if dims.averages && dims.subSpecs
            if sz(dims.subSpecs) >= 2
                ind = 1:size(MRS_struct.fids.data,2);
                ind = reshape(ind, [sz(dims.averages)*sz(dims.subSpecs)/sz(dims.subSpecs) sz(dims.subSpecs)])';
                ind = ind(:);
                MRS_struct.fids.data = MRS_struct.fids.data(:,ind);
            end
        end

    case 'Phillips'

end

end


function [dims, fids, npoints] = specify_dims(hdr, hdr_ext, fids)

% In NIfTI MRS, the three spatial dimensions and the time dimension occupy
% fixed indices in the (maximum) 7-D array
dims.x = 1;
dims.y = 2;
dims.z = 3;
dims.t = 4;

% There are some pre-defined dimension names according to the FID-A
% convention. These dimensions may or may not be stored in the NIfTI MRS
% header, so we'll initialize them as 0.
dims.coils    = 0;
dims.averages = 0;
dims.subSpecs = 0;
dims.extras   = 0;

% The NIfTI MRS standard reserves the remaining 3 dimensions, which are
% then explicitly specified in the JSON header extension fields dim_5,
% dim_6 and dim_7.
dims = parse_hdr_ext(hdr_ext, dims);

% Parse the NIfTI hdr.dim field:
all_dims = hdr.dim(2:end); % all dimensions (including singletons)

% Find the number of points
npoints = all_dims(dims.t);

% ORDERING THE DATA AND DIMENSIONS
% The FID-A array ordering conventions differ from the NIfTI MRS
% convention.
if prod(all_dims(1:3)) == 1 % x=y=z=1
    dims.x = 0;
    dims.y = 0;
    dims.z = 0;
    fids   = squeeze(fids);

    %Now that we've indexed the dimensions of the data array, we now need to
    %permute it so that the order of the dimensions is standardized:  we want
    %the order to be as follows:
    %   1) time domain data.
    %   2) coils.
    %   3) averages.
    %   4) subSpecs.
    %   5) extras.

    % Adjust dimension indices for the fact that we have collapsed the
    % three spatial dimensions (which we don't need for SVS data)
    sqz_dims = {};
    dims_fieldnames = fieldnames(dims);
    for jj = 1:length(dims_fieldnames)
        if dims.(dims_fieldnames{jj}) ~= 0
            % Subtract 3 (x, y, z) from the dimension indices
            dims.(dims_fieldnames{jj}) = dims.(dims_fieldnames{jj}) - 3;
            sqz_dims{end+1} = dims_fieldnames{jj}; %#ok<*AGROW>
        end
    end

    if length(sqz_dims) == 5
        fids          = permute(fids, [dims.t dims.coils dims.averages dims.subSpecs dims.extras]);
        dims.t        = 1;
        dims.coils    = 2;
        dims.averages = 3;
        dims.subSpecs = 4;
        dims.extras   = 5;
    elseif length(sqz_dims) == 4
        if dims.extras == 0
            fids          = permute(fids, [dims.t dims.coils dims.averages dims.subSpecs]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 3;
            dims.subSpecs = 4;
            dims.extras   = 0;
        elseif dims.subSpecs == 0
            fids          = permute(fids, [dims.t dims.coils dims.averages dims.extras]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 3;
            dims.subSpecs = 0;
            dims.extras   = 4;
        elseif dims.averages == 0
            fids          = permute(fids, [dims.t dims.coils dims.subSpecs dims.extras]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 0;
            dims.subSpecs = 3;
            dims.extras   = 4;
        elseif dims.coils == 0
            fids          = permute(fids, [dims.t dims.averages dims.subSpecs dims.extras]);
            dims.t        = 1;
            dims.coils    = 0;
            dims.averages = 2;
            dims.subSpecs = 3;
            dims.extras   = 4;
        end
    elseif length(sqz_dims) == 3
        if dims.extras == 0 && dims.subSpecs == 0
            fids = permute(fids, [dims.t dims.coils dims.averages]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 3;
            dims.subSpecs = 0;
            dims.extras   = 0;
        elseif dims.extras == 0 && dims.averages == 0
            fids          = permute(fids, [dims.t dims.coils dims.subSpecs]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 0;
            dims.subSpecs = 3;
            dims.extras   = 0;
        elseif dims.extras == 0 && dims.coils == 0
            fids          = permute(fids, [dims.t dims.averages dims.subSpecs]);
            dims.t        = 1;
            dims.coils    = 0;
            dims.averages = 2;
            dims.subSpecs = 3;
            dims.extras   = 0;
        end
    elseif length(sqz_dims) == 2
        if dims.extras == 0 && dims.subSpecs == 0 && dims.averages == 0
            fids          = permute(fids, [dims.t dims.coils]);
            dims.t        = 1;
            dims.coils    = 2;
            dims.averages = 0;
            dims.subSpecs = 0;
            dims.extras   = 0;
        elseif dims.extras == 0 && dims.subSpecs == 0 && dims.coils == 0
            fids          = permute(fids, [dims.t dims.averages]);
            dims.t        = 1;
            dims.coils    = 0;
            dims.averages = 2;
            dims.subSpecs = 0;
            dims.extras   = 0;
        elseif dims.extras == 0 && dims.averages == 0 && dims.coils == 0
            fids          = permute(fids, [dims.t dims.subSpecs]);
            dims.t        = 1;
            dims.coils    = 0;
            dims.averages = 0;
            dims.subSpecs = 2;
            dims.extras   = 0;
        end
    elseif length(sqz_dims) == 1
        dims.t        = 1;
        dims.coils    = 0;
        dims.averages = 0;
        dims.subSpecs = 0;
        dims.extras   = 0;
    end

else

    error(sprintf('\nData are not single-voxel data. Exiting...\n\n'));

end

end


function dims = parse_hdr_ext(hdr_ext, dims)

if isfield(hdr_ext, 'dim_5')
    dim_number = 5;
    % This field may come in as a cell or a string.
    if iscell(hdr_ext.dim_5)
        dim_5 = hdr_ext.dim_5{1};
    else
        dim_5 = hdr_ext.dim_5;
    end
    switch dim_5
        case 'DIM_COIL'
            dims.coils      = dim_number;
        case 'DIM_DYN'
            dims.averages   = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras     = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras     = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras     = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs   = dim_number;
        case 'DIM_MEAS'
            dims.extras     = dim_number;
        case 'DIM_USER_0'
            dims.extras     = dim_number;
        case 'DIM_USER_1'
            dims.extras     = dim_number;
        case 'DIM_USER_2'
            dims.extras     = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs   = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_5: %s', dim_5);
    end
end

if isfield(hdr_ext, 'dim_6')
    dim_number = 6;
    % This field may come in as a cell or a string.
    if iscell(hdr_ext.dim_6)
        dim_6 = hdr_ext.dim_6{1};
    else
        dim_6 = hdr_ext.dim_6;
    end
    switch dim_6
        case 'DIM_COIL'
            dims.coils    = dim_number;
        case 'DIM_DYN'
            dims.averages = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras   = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras   = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras   = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras   = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs = dim_number;
        case 'DIM_MEAS'
            dims.extras   = dim_number;
        case 'DIM_USER_0'
            dims.extras   = dim_number;
        case 'DIM_USER_1'
            dims.extras   = dim_number;
        case 'DIM_USER_2'
            dims.extras   = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_6: %s', dim_6);
    end
end

if isfield(hdr_ext, 'dim_7')
    dim_number = 7;
    % This field may come in as a cell or a string.
    if iscell(hdr_ext.dim_7)
        dim_7 = hdr_ext.dim_7{1};
    else
        dim_7 = hdr_ext.dim_7;
    end
    switch dim_7
        case 'DIM_COIL'
            dims.coils    = dim_number;
        case 'DIM_DYN'
            dims.averages = dim_number;
        case 'DIM_INDIRECT_0'
            dims.extras   = dim_number;
        case 'DIM_INDIRECT_1'
            dims.extras   = dim_number;
        case 'DIM_INDIRECT_2'
            dims.extras   = dim_number;
        case 'DIM_PHASE_CYCLE'
            dims.extras   = dim_number;
        case 'DIM_EDIT'
            dims.subSpecs = dim_number;
        case 'DIM_MEAS'
            dims.extras   = dim_number;
        case 'DIM_USER_0'
            dims.extras   = dim_number;
        case 'DIM_USER_1'
            dims.extras   = dim_number;
        case 'DIM_USER_2'
            dims.extras   = dim_number;
        case 'DIM_ISIS'
            dims.subSpecs = dim_number;
        otherwise
            error('Unknown dimension value specified in dim_7: %s', dim_7);
    end
end

end



