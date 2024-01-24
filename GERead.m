function MRS_struct = GERead(MRS_struct, fname)
% 141007: RTN edits to accommodate Noeske version
% 160916: MM & RTN edits to accommodate different encoding schemes
% 180404: RTN edits for more flexible handling of different P-file header
%         revisions; added support for rdbm_rev_num 26.002
% 201023: MM added support for rdbm_rev_num 27.x
% 230728: MM added support for rdbm_rev_num 30
% 231212: MM added support for receiver phase toggle for sLASER WIP (thanks
%         RTN)

ii = MRS_struct.ii;

fid = fopen(fname, 'r', 'ieee-be');
if fid == -1
    error('File ''%s'' not found!', fname);
end
fseek(fid, 0, 'bof');
rdbm_rev_num = fread(fid, 1, 'real*4');
if rdbm_rev_num == 7.0
    pfile_header_size = 39984; % LX
elseif rdbm_rev_num == 8.0
    pfile_header_size = 60464; % Cardiac / MGD
elseif rdbm_rev_num > 5.0 && rdbm_rev_num < 6.0
    pfile_header_size = 39940; % Signa 5.5
else
    % In 11.0 and later the header and data are stored as little-endian
    fclose(fid);
    fid = fopen(fname, 'r', 'ieee-le');
    fseek(fid, 0, 'bof');
    rdbm_rev_num = fread(fid, 1, 'real*4');
    if rdbm_rev_num == 9.0 % 11.0 product release
        pfile_header_size = 61464;
    elseif rdbm_rev_num == 11.0 % 12.0 product release
        pfile_header_size = 66072;
    end
end

MRS_struct.p.GE.rdbm_rev_num(ii) = rdbm_rev_num;
chkRev = {'14.3','16','20.006','20.007','24','26.002','27','27.001','28.002','28.003','30'};
assert(any(strcmp(num2str(rdbm_rev_num), chkRev)), ...
    sprintf(['GERead.m is not fully functional with P-file header revision number %g. ' ...
             'Please contact the Gannet developers (mam4041@med.cornell.edu) for assistance.'], ...
    rdbm_rev_num));

% RTN 2018
% Added flexible P-file revision support
% Values are read from rdb_hdr and image sub-headers
% Position can be found in rdbm.h (RDB_HEADER_REC) and imagedb.h (MRIMAGEDATATYPE)

% RTN 2018
% unsigned int rdb_hdr_ps_mps_freq
% float rdb_hdr_user0
% float rdb_hdr_user4
% float rdb_hdr_user19
% short rdb_hdr_nechoes
% short rdb_hdr_navs
% short rdb_hdr_nframes
% short rdb_hdr_point_size
% unsigned short rdb_hdr_da_xres
% short rdb_hdr_da_yres
% short rdb_hdr_dab[0].start_rcv
% short rdb_hdr_dab[0].stop_rcv
% int rdb_hdr_off_image
% int rdb_hdr_off_data
%
% image sub-header
% int te
% int tr
% float user8-10    voxel dimensions
% float user19      rf waveform
% float user20-21   offset frequencies
% float user22      pulse width (-1 default)

switch num2str(rdbm_rev_num)

    case '14.3'

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 181;
        image_tr = 179;

        % float
        image_user8  = 38;
        image_user19 = 49;
        image_user20 = 50;
        image_user22 = 52;
        image_user24 = 56;

    case '16'

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 193;
        image_tr = 191;

        % float
        image_user8  = 50;
        image_user19 = 61;
        image_user20 = 62;
        image_user22 = 64;
        image_user24 = 68;

    case {'20.006','20.007','24'}

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 267;
        image_tr = 265;

        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;
        image_user24 = 116;

    case {'26.002','27','27.001','28.002','28.003','30'}

        % int
        rdb_hdr_off_image   = 11;
        rdb_hdr_off_data    = 2;
        rdb_hdr_ps_mps_freq = 123;

        % float
        rdb_hdr_user0  = 71;
        rdb_hdr_user4  = 75;
        rdb_hdr_user19 = 90;

        % short
        rdb_hdr_nechoes       = 74;
        rdb_hdr_navs          = 75;
        rdb_hdr_nframes       = 76;
        rdb_hdr_point_size    = 80;
        rdb_hdr_da_xres       = 90;
        rdb_hdr_da_yres       = 91;
        rdb_hdr_dab_start_rcv = 133;
        rdb_hdr_dab_stop_rcv  = 134;

        % int
        image_te = 267;
        image_tr = 265;

        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;
        image_user24 = 116;

end

% Read rdb header as short, int and float
fseek(fid, 0, 'bof');
hdr_value = fread(fid, rdb_hdr_dab_stop_rcv, 'integer*2');
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, rdb_hdr_user19, 'real*4');
fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');

if rdbm_rev_num > 11.0
    pfile_header_size = i_hdr_value(rdb_hdr_off_data);
end

MRS_struct.p.LarmorFreq(ii) = i_hdr_value(rdb_hdr_ps_mps_freq)/1e7;
MRS_struct.p.sw(ii)         = f_hdr_value(rdb_hdr_user0);
nechoes                     = hdr_value(rdb_hdr_nechoes);
MRS_struct.p.GE.nechoes(ii) = nechoes;
nex                         = hdr_value(rdb_hdr_navs);
MRS_struct.p.GE.NEX(ii)     = nex;
nframes                     = hdr_value(rdb_hdr_nframes);
point_size                  = hdr_value(rdb_hdr_point_size);
MRS_struct.p.npoints(ii)    = hdr_value(rdb_hdr_da_xres);
MRS_struct.p.nrows(ii)      = hdr_value(rdb_hdr_da_yres);

start_recv = hdr_value(rdb_hdr_dab_start_rcv);
stop_recv  = hdr_value(rdb_hdr_dab_stop_rcv);
nreceivers = (stop_recv - start_recv) + 1;

% RTN 2018
dataframes = f_hdr_value(rdb_hdr_user4)/nex;
refframes  = f_hdr_value(rdb_hdr_user19);

% Read image header as int and float
% Find TE/TR
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
t_hdr_value         = fread(fid, image_te, 'integer*4');
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
o_hdr_value         = fread(fid, image_user24, 'real*4');
MRS_struct.p.TE(ii) = t_hdr_value(image_te)/1e3;
MRS_struct.p.TR(ii) = t_hdr_value(image_tr)/1e3;

% Find voxel dimensions and edit pulse parameters
MRS_struct.p.voxdim(ii,:)             = o_hdr_value(image_user8:image_user8+2)';
MRS_struct.p.GE.editRF.waveform(ii)   = o_hdr_value(image_user19);
MRS_struct.p.GE.editRF.freq_Hz(ii,:)  = o_hdr_value(image_user20:image_user20+1)';
MRS_struct.p.GE.editRF.freq_ppm(ii,:) = (MRS_struct.p.GE.editRF.freq_Hz(ii,:) / MRS_struct.p.LarmorFreq(ii)) + 4.68;
MRS_struct.p.GE.editRF.dur(ii)        = o_hdr_value(image_user22)/1e3;
% RTN 2018: check for default value (-1) of pulse length
if MRS_struct.p.GE.editRF.dur(ii) <= 0
    MRS_struct.p.GE.editRF.dur(ii) = 16;
end
% CV24: this CV is especially important for the WIP HERMES/sLASER implementations
MRS_struct.p.GE.cv24(ii) = o_hdr_value(image_user24);

% Spectro prescan pfiles
if MRS_struct.p.npoints(ii) == 1 && MRS_struct.p.nrows(ii) == 1
    MRS_struct.p.npoints(ii) = 2048;
end

% Compute size (in bytes) of data
data_elements          = MRS_struct.p.npoints(ii) * 2;
totalframes            = MRS_struct.p.nrows(ii) * nechoes; % RTN nechoes mulitply
MRS_struct.p.nrows(ii) = totalframes;
data_elements          = data_elements * totalframes * nreceivers;

fseek(fid, pfile_header_size, 'bof');
% Read data: point_size = 2 means 16-bit data, point_size = 4 means EDR
if point_size == 2
    raw_data = fread(fid, data_elements, 'integer*2');
else
    raw_data = fread(fid, data_elements, 'integer*4');
end
fclose(fid);

% 110303 CJE
% Calculate Navg from nframes, 8 water frames, 2 phase cycles
% Needs to be specific to single experiment - for frame rejection
% RTN edits to accommodate Noeske version raee 20141007
% MM (160916): Incorporating more edits from RTN to handle dual-echo data
%              acquired with one of four possible encoding schemes:
%              NEX=2/noadd=0, NEX=2/noadd=1, NEX=8/noadd=0, NEX=8/noadd=1
% MM (171120): RTN edits to accomodate HERMES aquisitions; better looping
%              over phase cycles
% MM (200713): RTN edits for better handling of data if nechoes == 1
if nechoes == 1

    if (dataframes + refframes) ~= nframes
        mult                      = 1;
        MRS_struct.p.GE.noadd(ii) = 1;
        dataframes                = dataframes * nex;
        refframes                 = nframes - dataframes;
    else
        mult                      = 1/nex;
        MRS_struct.p.GE.noadd(ii) = 0;
    end

    ShapeData = reshape(raw_data, [2 MRS_struct.p.npoints(ii) totalframes nreceivers]);
    WaterData = ShapeData(:,:,2:refframes+1,:) * mult;
    MetabData = ShapeData(:,:,refframes+2:end,:) * mult / 2;

    totalframes = totalframes - (refframes + 1);

    MRS_struct.p.nrows(ii)       = totalframes;
    MRS_struct.p.nrows_water(ii) = refframes;
    MRS_struct.p.Navg(ii)        = dataframes * nex;
    MRS_struct.p.Nwateravg(ii)   = refframes * nex;

else

    MRS_struct.p.Navg(ii) = dataframes * nex * nechoes; % RTN 2017

    if (dataframes + refframes) ~= nframes
        mult                      = nex/2; % RTN 2016   1; % RTN 2017
        multw                     = nex;   % RTN 2016   1; % RTN 2017
        MRS_struct.p.GE.noadd(ii) = 1;
        dataframes                = dataframes * nex;
        refframes                 = nframes - dataframes;
    else
        mult                      = 1/nex; % MM 2020   1; % RTN 2017       nex/2; % RTN 2016
        multw                     = 1;     % MM 2020   1/nex; % RTN 2017   1; % RTN 2016
        MRS_struct.p.GE.noadd(ii) = 0;
    end

    MRS_struct.p.Nwateravg(ii)   = refframes * nechoes; % RTN 2017
    MRS_struct.p.nrows_water(ii) = refframes;

    if totalframes ~= (dataframes + refframes + 1) * nechoes % RTN 2017
        error('Number of totalframes does not equal (dataframes + refframes + 1) * nechoes');
    end

    ShapeData = reshape(raw_data, [2 MRS_struct.p.npoints(ii) totalframes nreceivers]);

    [X1,X2]   = ndgrid(1:refframes, 1:nechoes);
    X1        = X1'; X1 = X1(:);
    X2        = X2'; X2 = X2(:);
    Y1        = (-1).^(MRS_struct.p.GE.noadd(ii) * (X1-1));
    if MRS_struct.p.GE.cv24(ii) >= 16384 % Do not apply any phase cycling correction when the receiver phase toggle in sLASER has been set
        Y1    = ones(size(Y1,1),1);
    end
    Y1        = permute(repmat(Y1, [1 MRS_struct.p.npoints(ii) 2 nreceivers]), [3 2 1 4]);
    Y2        = 1 + (totalframes/nechoes) * (X2-1) + X1;
    WaterData = Y1 .* ShapeData(:,:,Y2,:) * multw;

    [X1,X2]   = ndgrid(1:dataframes, 1:nechoes);
    X1        = X1'; X1 = X1(:);
    X2        = X2'; X2 = X2(:);
    Y1        = (-1).^(MRS_struct.p.GE.noadd(ii) * (X1-1));
    if MRS_struct.p.GE.cv24(ii) >= 16384 % Do not apply any phase cycling correction when the receiver phase toggle in sLASER has been set
        Y1    = ones(size(Y1,1),1);
    end
    Y1        = permute(repmat(Y1, [1 MRS_struct.p.npoints(ii) 2 nreceivers]), [3 2 1 4]);
    Y2        = 1 + refframes + (totalframes/nechoes) * (X2-1) + X1;
    MetabData = Y1 .* ShapeData(:,:,Y2,:) * mult;

    totalframes            = totalframes - (refframes + 1) * nechoes; % RTN 2017
    MRS_struct.p.nrows(ii) = totalframes;

end

MetabData = squeeze(complex(MetabData(1,:,:,:), MetabData(2,:,:,:)));
MetabData = permute(MetabData, [3 1 2]);
if size(MetabData,1) == 1 % re-permute array dimensions in cases were there is only one average
    MetabData = permute(MetabData, [3 2 1]);
end
WaterData = squeeze(complex(WaterData(1,:,:,:), WaterData(2,:,:,:)));
WaterData = permute(WaterData, [3 1 2]);
if size(WaterData,1) == 1 % re-permute array dimensions in cases were there is only one average
    WaterData = permute(WaterData, [3 2 1]);
end

% Combine coils using generalized least squares method (An et al., JMRI,
% 2013, doi:10.1002/jmri.23941); the noise covariance matrix is more
% optimally estimated by using all averages as suggested by Rodgers &
% Robson (MRM, 2010, doi:10.1002/mrm.22230)
[nCh, nPts, nReps]             = size(WaterData);
noise_pts                      = false(1,nPts);
noise_pts(ceil(0.75*nPts):end) = true;
noise_pts                      = repmat(noise_pts, [1 nReps]);
tmpWaterData                   = reshape(WaterData, [nCh nPts*nReps]);

e             = tmpWaterData(:,noise_pts);
Psi           = e*e';
WaterData_avg = mean(WaterData,3);
S             = WaterData_avg(:,1);
w             = (S'*(Psi\S))^-1 * S' / Psi;
WaterData     = w.' .* WaterData;

MRS_struct.fids.data_water = mean(squeeze(sum(WaterData,1)),2);

[nCh, nPts, nReps]             = size(MetabData);
noise_pts                      = false(1,nPts);
noise_pts(ceil(0.75*nPts):end) = true;
noise_pts                      = repmat(noise_pts, [1 nReps]);
tmpMetabData                   = reshape(MetabData, [nCh nPts*nReps]);

e         = tmpMetabData(:,noise_pts);
Psi       = e*e';
w         = (S'*(Psi\S))^-1 * S' / Psi;
MetabData = w.' .* MetabData;

MRS_struct.fids.data = squeeze(sum(MetabData,1));

end



