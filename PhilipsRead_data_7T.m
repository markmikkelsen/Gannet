function [ MRS_struct ] = PhilipsRead_data_7T(MRS_struct, fname, fname_water )
% RE/CJE Parse SPAR file for header info
% 110825

%ADH 20130701 - modify such that if the water reference scan is loaded in
%the same file, it works. for this case, will only have one input fname
%because all the raw data is in the same place. Next step will be to put
%a flag to see if there is water data in a separate scan, in the same scan
%as a water ref (as called in SDAT) or if its just not been collected at all

% will need to get all the indicies first because will need to know how big
% the waterdata is before loading in the GABAdata.

%% Read in water data --  MGSaleh 2017
if nargin > 2 
 % work out data header name
sparname = [fname_water(1:(end-4)) 'list'];
sparheader = textscan(sparname, '%s');
%ADH - decide if there is water data as ref data included in the data and
%if so, set a flag to pull it out properly...
sparidx=find(ismember(sparheader, 'F-resolution')==1);
MRS_struct.p.npoints_water = str2double(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1);
MRS_struct.p.nrows_water = str2double(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1);
MRS_struct.p.NSA_water(MRS_struct.ii) = str2double(sparheader{sparidx+2});
MRS_struct.p.Navg_water(MRS_struct.ii) = MRS_struct.p.NSA_water(MRS_struct.ii)*MRS_struct.p.nrows_water;
%Need to determine the offset of the data - i.e. how many channels are
%there...
sparidx=find(ismember(sparheader, 'NOI')==1);
MRS_struct.p.coil_channels=size(sparidx,1)-2;
sparidx=find(ismember(sparheader, 'STD')==1);
MRS_struct.p.ptr_offset_water=str2double(sparheader{sparidx(3)+20});

MRS_struct.p.Navg_water(MRS_struct.ii) = MRS_struct.p.Navg_water(MRS_struct.ii)*MRS_struct.p.coil_channels; % MGSaleh
%Need to skip rows associated with the '
%                                                        [ %real/imag FID points random total_FIDS/dynamic scans dynamic scans]
MRS_struct.fids.data_water = readraw_Gannet(fname_water, 'float', [2 MRS_struct.p.npoints_water 1 MRS_struct.p.Navg_water(MRS_struct.ii)/MRS_struct.p.nrows_water MRS_struct.p.nrows_water], 'l',MRS_struct.p.ptr_offset_water);
%  Make data complex.
MRS_struct.fids.data_water = squeeze(complex(MRS_struct.fids.data_water(1,:,:,:,:), MRS_struct.fids.data_water(2,:,:,:,:)));
FullData_water = MRS_struct.fids.data_water;

% Coil combination -- MM and MGSaleh 2017
firstpoint = mean(conj(FullData_water(1,:,:)),3);
% firstpoint_water = firstpoint;
channels_scale = squeeze(sqrt(sum(firstpoint.*conj(firstpoint))));
firstpoint = repmat(firstpoint,[MRS_struct.p.npoints_water 1 size(FullData_water,3)])/channels_scale;
FullData_water = FullData_water .* firstpoint;  %Multiply FIDs to the first point of the FID -- MGsaleh 2017
FullData_water = conj(squeeze(sum(FullData_water,2))); %Combine all channels by simply adding them -- MGSaleh 2017
MRS_struct.fids.data_water=FullData_water;
%undo time series phase cycling to match GE
corrph = ones(size(MRS_struct.fids.data_water));
for jj=1:size(MRS_struct.fids.data_water,2)
    corrph(:,jj) = corrph(:,jj) * (-1).^(jj+1);
end
MRS_struct.fids.data_water = MRS_struct.fids.data_water .* corrph;  %Flips the real-time data about x-axis -- MGSaleh 2017
%Philips data appear to be phased already (ideal case) %Code below phases spectra to ensure everything is upright -- MGSaleh 2017
MRS_struct.fids.data_water = MRS_struct.fids.data_water .*repmat(conj(MRS_struct.fids.data_water(1,:))./abs(MRS_struct.fids.data_water(1,:)),[MRS_struct.p.npoints_water 1]); %Phases spectra to ensure everything is upright -- MGSaleh 2017
disp('water data ... done')
end
%% Read in water suppressed data --  MGSaleh 2017
% work out data header name
sparname = [fname(1:(end-4)) 'list'];
sparheader = textscan(sparname, '%s');
%ADH - decide if there is water data as ref data included in the data and
%if so, set a flag to pull it out properly...
% sparidx=find(ismember(sparheader, 'number_of_mixes')==1);
% nu_mixes = str2double(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'F-resolution')==1);
MRS_struct.p.npoints = str2double(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'number_of_extra_attribute_1_values')==1); %Dynamic scans
MRS_struct.p.nrows = str2double(sparheader{sparidx+2});
sparidx=find(ismember(sparheader, 'number_of_signal_averages')==1);
MRS_struct.p.Navg(MRS_struct.ii) = str2double(sparheader{sparidx+2});
%Need to determine the offset of the data - i.e. how many channels are there...
sparidx=find(ismember(sparheader, 'NOI')==1);
MRS_struct.p.coil_channels=size(sparidx,1)-2;
sparidx=find(ismember(sparheader, 'STD')==1);
MRS_struct.p.ptr_offset=str2double(sparheader{sparidx(3)+20});
MRS_struct.p.Navg_all_chann(MRS_struct.ii) = MRS_struct.p.Navg(MRS_struct.ii)*MRS_struct.p.nrows*MRS_struct.p.coil_channels; % MGSaleh
%Need to skip rows associated with the                [real/imag   FID points             random      total_FIDS/dynamic_scans                                            dynamic scans  ]
MRS_struct.fids.data = readraw_Gannet(fname, 'float', [    2       MRS_struct.p.npoints      1        MRS_struct.p.Navg_all_chann(MRS_struct.ii)/MRS_struct.p.nrows    MRS_struct.p.nrows], 'l',MRS_struct.p.ptr_offset);
%  Make data complex.
MRS_struct.fids.data = squeeze(complex(MRS_struct.fids.data(1,:,:,:,:), MRS_struct.fids.data(2,:,:,:,:)));
%                                                   [     FID points             coil          NSA            dynamic scans   ]
MRS_struct.fids.data = reshape(MRS_struct.fids.data,[size(MRS_struct.fids.data,1) 32   MRS_struct.p.Navg    MRS_struct.p.nrows]);
MRS_struct.fids.data = reshape(MRS_struct.fids.data, [size(MRS_struct.fids.data,1) size(MRS_struct.fids.data,2) ...
                       size(MRS_struct.fids.data,3)*size(MRS_struct.fids.data,4)]); % Combining NSA and dynamic scans -- MGSaleh 2017
FullData = MRS_struct.fids.data;

% Coil combination -- MM and MGSaleh 2017
firstpoint = mean(conj(FullData(1,:,:)),3);
channels_scale = squeeze(sqrt(sum(firstpoint.*conj(firstpoint))));
firstpoint = repmat(firstpoint,[MRS_struct.p.npoints 1 size(FullData,3)])/channels_scale;
FullData = FullData .* firstpoint;  %Multiply FIDs to the first point of the FID -- MGsaleh 2017
FullData = conj(squeeze(sum(FullData,2))); %Combine all channels by simply adding them -- MGSaleh 2017
MRS_struct.fids.data=FullData;
corrph = repmat([ones(1,MRS_struct.p.Navg) -1*ones(1,MRS_struct.p.Navg)], [1 size(MRS_struct.fids.data,2)/(MRS_struct.p.Navg*2)]);
corrph = repmat(corrph, [size(MRS_struct.fids.data,1) 1]);
MRS_struct.fids.data = MRS_struct.fids.data .* corrph;  %Flips the real-time data about x-axis -- MGSaleh 2017
if MRS_struct.p.HERMES % MM (170604): Use NAA signal for zero-order pre-phasing
    inds = [ones(1,4), zeros(1,4), zeros(1,4), ones(1,4)];
    inds = repmat(inds, [1 size(MRS_struct.fids.data,2)/numel(inds)]);
    phi = conj(MRS_struct.fids.data(1,inds==1))./abs(MRS_struct.fids.data(1,inds==1));
    phi = repelem(phi,2);
    MRS_struct.fids.data = MRS_struct.fids.data .* repmat(phi,[MRS_struct.p.npoints 1]);
else
    MRS_struct.fids.data = MRS_struct.fids.data .* repmat(conj(MRS_struct.fids.data(1,:))./abs(MRS_struct.fids.data(1,:)),[MRS_struct.p.npoints 1]);
end
disp('water suppressed data ... done')
%I moved it to the end of the function -- MGSaleh 05252018
MRS_struct.p.Navg(MRS_struct.ii) = MRS_struct.p.Navg(MRS_struct.ii)*MRS_struct.p.nrows;

end

