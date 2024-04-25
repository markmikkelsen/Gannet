function MRS_struct = PhilipsRead(MRS_struct, fname, fname_water)
% RE/CJE Parse SPAR file for header info

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PhilipsRead is designed to handle 'off-first' data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = MRS_struct.ii;

% Work out data header name
[~,~,ext] = fileparts(fname);
if all(isstrprop(ext(end-3:end), 'upper'))
    spar_ext = 'SPAR';
else
    spar_ext = 'spar';
end
sparname = [fname(1:(end-4)) spar_ext];
sparname = fopen(sparname,'r');
sparheader = textscan(sparname, '%s');
sparheader = sparheader{1};
sparidx = find(ismember(sparheader, 'samples') == 1);
MRS_struct.p.npoints(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'rows') == 1);
MRS_struct.p.nrows(ii) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'averages') == 1);
if strcmpi(MRS_struct.p.seqorig, 'Philips')
    MRS_struct.p.Navg(ii) = MRS_struct.p.nrows(ii);
else
    MRS_struct.p.Navg(ii) = MRS_struct.p.nrows(ii) * str2double(sparheader{sparidx+2});
end
sparidx = find(ismember(sparheader, 'repetition_time') == 1);
MRS_struct.p.TR(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'echo_time') == 1);
MRS_struct.p.TE(ii) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'synthesizer_frequency') == 1);
MRS_struct.p.LarmorFreq(ii) = str2double(sparheader{sparidx+2})/1e6;
sparidx = find(ismember(sparheader, 'sample_frequency') == 1);
MRS_struct.p.sw(ii) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_size') == 1);
MRS_struct.p.voxdim(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_size') == 1);
MRS_struct.p.voxdim(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_size') == 1);
MRS_struct.p.voxdim(ii,3) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_off_center') == 1);
MRS_struct.p.voxoff(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_off_center') == 1);
MRS_struct.p.voxoff(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_off_center') == 1);
MRS_struct.p.voxoff(ii,3) = str2double(sparheader{sparidx+2});

sparidx = find(ismember(sparheader, 'ap_angulation') == 1);
MRS_struct.p.voxang(ii,2) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'lr_angulation') == 1);
MRS_struct.p.voxang(ii,1) = str2double(sparheader{sparidx+2});
sparidx = find(ismember(sparheader, 'cc_angulation') == 1);
MRS_struct.p.voxang(ii,3) = str2double(sparheader{sparidx+2});

MRS_struct.fids.data = SDATreadMEGA(fname, MRS_struct.p.npoints(ii), MRS_struct.p.nrows(ii));

% Undo phase cycling
corrph = repmat([-1 1], [1 size(MRS_struct.fids.data,2)/2]);
corrph = repmat(corrph, [size(MRS_struct.fids.data,1) 1]);
MRS_struct.fids.data = MRS_struct.fids.data .* corrph;

% Re-introduce initial phase step
if MRS_struct.p.HERMES
%     if strcmp(MRS_struct.p.ON_OFF_order,'offfirst')
        phi       = repelem(conj(MRS_struct.fids.data(1,2:2:end)) ./ abs(MRS_struct.fids.data(1,2:2:end)),2);
%     elseif strcmp(MRS_struct.p.ON_OFF_order,'onfirst')
%         ind1      = sort([1:4:size(MRS_struct.fids.data,2) 2:4:size(MRS_struct.fids.data,2)]);
%         ind2      = sort([3:4:size(MRS_struct.fids.data,2) 4:4:size(MRS_struct.fids.data,2)]);
%         phi(ind1) = repelem(conj(MRS_struct.fids.data(1,1:4:end)) ./ abs(MRS_struct.fids.data(1,1:4:end)),2);
%         phi(ind2) = repelem(conj(MRS_struct.fids.data(1,4:4:end)) ./ abs(MRS_struct.fids.data(1,4:4:end)),2);
%     end
    MRS_struct.fids.data = MRS_struct.fids.data .* repmat(phi, [MRS_struct.p.npoints(ii) 1]);
else
    if strcmp(MRS_struct.p.target{1}, 'GSH')
        MRS_struct.fids.data = MRS_struct.fids.data .* ...
            repmat(conj(mean(MRS_struct.fids.data(1,:))) ./ abs(mean(MRS_struct.fids.data(1,:))), size(MRS_struct.fids.data));
    else
        MRS_struct.fids.data = MRS_struct.fids.data .* ...
            repmat(conj(MRS_struct.fids.data(1,:)) ./ abs(MRS_struct.fids.data(1,:)), [MRS_struct.p.npoints(ii) 1]);
    end
end
% Philips data appear to be phased already (ideal case)

MRS_struct.fids.data = conj(MRS_struct.fids.data);

if nargin > 2
    % Load water data
    MRS_struct.p.Nwateravg       = 1; % water SDAT is average not sum
    MRS_struct.p.nrows_water     = 1;
    MRS_struct.fids.data_water   = SDATread(fname_water, MRS_struct.p.npoints(ii));
    MRS_struct.fids.data_water   = MRS_struct.fids.data_water .* ...
                                   conj(MRS_struct.fids.data_water(1)) ./ abs(MRS_struct.fids.data_water(1));
    MRS_struct.fids.data_water   = MRS_struct.fids.data_water.';
end

end



