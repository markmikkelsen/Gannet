function MRS_struct = GannetLoad(varargin)
% Gannet 3
% Created by RAEE (Nov. 5, 2012)
% Updates by MM, GO, MGS (2016-2023)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workflow summary
%   1. Pre-initialise
%   2. Determine data parameters from headers
%   3. Load data from files
%   4. Reconstruction of coil-sensitivity maps (PRIAM only)
%   5. Apply appropriate pre-processing
%   6. Build GannetLoad output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

MRS_struct.version.Gannet = '3.3.2';
MRS_struct.version.load   = '231212';
VersionCheck(0, MRS_struct.version.Gannet);
ToolboxCheck;

fprintf('\nGannet v%s - %s\n', MRS_struct.version.Gannet, ...
    hyperlink('https://github.com/markmikkelsen/Gannet', ...
    'https://github.com/markmikkelsen/Gannet'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -1. Check if GUI version of Gannet was run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gui_flag = false;
var_args = varargin;
num_args = nargin;

if islogical(var_args{1}) % first argument will be true if this is called from the GUI
    gui_flag    = true; % let the script know that the GUI was run, so we should look for a configuration file
    config_path = var_args{2}; % store the path to the configuration file (passed in as second argument)
    var_args    = var_args(3:end); % delete the guiFlag and config_path from var_args so as not to impact the rest of the script
    num_args    = length(var_args); % change nargin to account for the fact that we just removed the two GUI-related variables from var_args
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0. Parse the input arguments and check for typos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metabfile = var_args{1};
metabfile = GetFullPath(metabfile);
missing = 0;
for filecheck = 1:numel(metabfile)
    if ~exist(metabfile{filecheck}, 'file')
        fprintf('\nThe file ''%s'' (#%d) is missing. Typo?\n', metabfile{filecheck}, filecheck);
        missing = 1;
    end
end

if num_args > 1 && ~isempty(var_args{2})
    waterfile = var_args{2};
    waterfile = GetFullPath(waterfile);
    for filecheck = 1:numel(waterfile)
        if ~exist(waterfile{filecheck}, 'file')
            fprintf('\nThe water reference file ''%s'' (#%d) is missing. Typo?\n', waterfile{filecheck}, filecheck);
            missing = 1;
        end
    end
    [~,~,ext] = fileparts(metabfile{1});
    if ~strcmpi(ext, '.rda')
        assert(isequal(size(metabfile), size(waterfile)), 'The metabolite and water reference filename cell array inputs must have the same dimensions.');
    end
end

if missing
    fprintf('\n');
    error('Not all input files could be found. Please check filenames. Exiting...');
end

if num_args == 3
    MRS_struct.p.trim_avgs = 1;
    if isnumeric(var_args{3})
        trimAvgs = var_args{3};
        assert(size(trimAvgs,2) == 2, 'The third input argument must be a M x 2 array.');
    else
        [~,~,ext] = fileparts(var_args{3});
        assert(strcmpi(ext, '.csv'), 'The third argument must be a .csv file.');
        trimAvgs = readtable(var_args{3});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if gui_flag % if we launched this script from the GUI, have GannetPreInitialise read from a scan configuration file created by the GUI
    MRS_struct = GannetPreInitialiseGUIVersion(config_path, MRS_struct);
else % otherwise, run GannetPreInitialise as usual
    MRS_struct = GannetPreInitialise(MRS_struct);
end

CheckTargets(MRS_struct);

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end

if (size(metabfile,2) == 1 && MRS_struct.p.join == 0) || (MRS_struct.p.join == 1 && size(metabfile,1) == 1)
    metabfile = metabfile';
elseif MRS_struct.p.join == 0 && any(size(metabfile) > 1)
    metabfile = metabfile(:)';
end
MRS_struct.metabfile = metabfile;

if exist('waterfile', 'var')
    [~,~,ext] = fileparts(waterfile{1});
    if (size(waterfile,2) == 1 && MRS_struct.p.join == 0) || (MRS_struct.p.join == 1 && size(waterfile,1) == 1 && ~strcmpi(ext, '.rda'))
        waterfile = waterfile';
    elseif MRS_struct.p.join == 0 && any(size(waterfile) > 1)
        waterfile = waterfile(:)';
    end
    MRS_struct.waterfile = waterfile;
end

if MRS_struct.p.phantom
    if MRS_struct.p.HERMES
        out = input('What was the order of the HERMES editing pulses in the experiment? E.g., CBAD: ','s');
    else
        out = input('Which editing pulse was first in the experiment? ON or OFF: ','s');
    end
    if strcmpi(out, 'OFF')
        MRS_struct.p.ON_OFF_order = 'offfirst';
    elseif strcmpi(out, 'ON')
        MRS_struct.p.ON_OFF_order = 'onfirst';
    else
        MRS_struct.p.ON_OFF_order = out;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discern input data format
MRS_struct = DiscernDataType(metabfile{1}, MRS_struct);

% Determine number of provided water-suppressed files in the batch
[numFilesPerScan, numScans] = size(metabfile);

% For Siemens RDA, each acquisition has two files so correct the number
if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
    numScans = numScans/2;
end

if MRS_struct.p.join
    fprintf('\nRunning GannetLoad in ''join'' mode...\n');
end

% Determine number of provided water-unsuppressed files in the batch
if exist('waterfile', 'var')
    MRS_struct.p.reference = 'H2O';
    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
        assert(size(metabfile,2)/2 == size(waterfile,2), ...
            'Number of single water-unsuppressed RDA files per subject does not match number of paired (i.e., ON/OFF) water-suppressed RDA files per subject.');
    else
        assert(size(metabfile,2) == size(waterfile,2), 'Number of water-unsuppressed files per subject does not match number of water-suppressed files per subject.');
    end
else
    MRS_struct.p.reference = 'Cr';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.p.numScans        = numScans;
MRS_struct.p.numFilesPerScan = numFilesPerScan;
run_count    = 0;
error_report = cell(1);
catch_ind    = 1;

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','MATLAB:rankDeficientMatrix');

for ii = 1:MRS_struct.p.numScans % Loop over all files in the batch (from metabfile)
    
    [~,b,c] = fileparts(metabfile{1,ii});
    if MRS_struct.p.join
        fprintf('\nLoading %s and %i other file(s)...', [b c], numFilesPerScan - 1);
    else
        fprintf('\nLoading %s...', [b c]);
    end
    
    f = dir(metabfile{1,ii});
    if f.bytes/1e6 > 150
        fprintf('\nLarge file detected (%.2f MB). Please wait...', f.bytes/1e6);
    end
    
    try % pass to next dataset if errors occur
        
        MRS_struct.ii = ii;
        
        if strcmp(MRS_struct.p.vendor, 'Philips_raw')
            
            MRS_struct = PhilipsRawRead(MRS_struct, metabfile{1,ii}, 3, 0);
            MRS_struct.fids.data = conj(squeeze(MRS_struct.multivoxel.allsignals(:,:,1,:)));
            if exist('waterfile', 'var')
                MRS_struct.p.reference = 'H2O';
            end
            
        else
            
            switch MRS_struct.p.vendor
                case 'DICOM'
                    loadFun = @DICOMRead;
                case 'GE'
                    loadFun = @GERead;
                case 'NIfTI'
                    loadFun = @NIfTIMRSRead;
                case 'Philips'
                    loadFun = @PhilipsRead;
                case 'Philips_data'
                    loadFun = @PhilipsDataRead;
                case 'Siemens_dicom'
                    loadFun = @SiemensDICOMRead;
                case 'Siemens_rda'
                    loadFun = @SiemensRead;
                case 'Siemens_twix'
                    loadFun = @SiemensTWIXRead;
            end
            
            if exist('waterfile', 'var')
                if MRS_struct.p.join
                    MRS_struct = loadFun(MRS_struct, metabfile{1,ii}, waterfile{1,ii});
                    for kk = 2:numFilesPerScan
                        sub_MRS_struct             = loadFun(MRS_struct, metabfile{kk,ii}, waterfile{kk,ii});
                        MRS_struct.fids.data       = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                        MRS_struct.fids.data_water = [MRS_struct.fids.data_water sub_MRS_struct.fids.data_water];
                    end
                    MRS_struct.p.nrows(ii) = size(MRS_struct.fids.data,2);
                    MRS_struct.p.Navg(ii)  = size(MRS_struct.fids.data,2);
                else
                    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                        MRS_struct = loadFun(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1}, waterfile{ii});
                    else
                        MRS_struct = loadFun(MRS_struct, metabfile{ii}, waterfile{ii});
                    end
                end
            else
                if MRS_struct.p.join
                    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                        MRS_struct = loadFun(MRS_struct, metabfile{1,ii*2}, metabfile{1,ii*2-1});
                    else
                        MRS_struct = loadFun(MRS_struct, metabfile{1,ii});
                    end
                    for kk = 2:numFilesPerScan
                        if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                            sub_MRS_struct = loadFun(MRS_struct, metabfile{kk,ii*2}, metabfile{kk,ii*2-1});
                        else
                            sub_MRS_struct = loadFun(MRS_struct, metabfile{kk,ii});
                        end
                        MRS_struct.fids.data = [MRS_struct.fids.data sub_MRS_struct.fids.data];
                    end
                    MRS_struct.p.nrows(ii) = size(MRS_struct.fids.data,2);
                    MRS_struct.p.Navg(ii)  = size(MRS_struct.fids.data,2);
                else
                    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                        MRS_struct = loadFun(MRS_struct, metabfile{ii*2}, metabfile{ii*2-1});
                    else
                        MRS_struct = loadFun(MRS_struct, metabfile{ii});
                    end
                end
            end
            
        end
        
        switch MRS_struct.p.vendor
            case 'GE'
                MRS_struct.p.reference = 'H2O';
            case 'Philips_data'
                if isfield(MRS_struct.fids, 'data_water')
                    MRS_struct.p.reference = 'H2O';
                else
                    MRS_struct.p.reference = 'Cr';
                end
        end
        
        % If user wants to trim scans
        if isfield(MRS_struct.p, 'trim_avgs')
            if isnumeric(var_args{3})
                if size(trimAvgs,1) > 1
                    t_start = trimAvgs(ii,1);
                    t_end   = trimAvgs(ii,2);
                else
                    t_start = trimAvgs(1,1);
                    t_end   = trimAvgs(1,2);
                end
            else
                fileInd = find(strcmpi(metabfile(ii), trimAvgs{:,1}));
                t_start = trimAvgs{fileInd,2};
                t_end   = trimAvgs{fileInd,3};
            end
            % Check if t_end is less than or equal to the total number of
            % acquired averages
            assert(t_end <= size(MRS_struct.fids.data,2), 'The requested trim exceeds the total number of acquired averages.');
            % Make sure t_start is an odd number and t_end is an even number
            if ~rem(t_start,2)
                t_start = t_start - 1;
            end
            if rem(t_end,2)
                t_end = t_end + 1;
            end
            MRS_struct.fids.data  = MRS_struct.fids.data(:,t_start:t_end);
            MRS_struct.p.Navg(ii) = size(MRS_struct.fids.data,2);
        end
        
        % Determine ON/OFF order
        MRS_struct = SpecifyOnOffOrder(MRS_struct);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   4. Reconstruction of coil-sensitivity maps
        %      (PRIAM only)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % if a PRIAM dataset is processed, load the coil reference scan and
        % calculate the SENSE reconstruction matrix here
        if MRS_struct.p.PRIAM
            MRS_struct = senseRecon(MRS_struct);
            PRIAMData = zeros(length(MRS_struct.p.vox), MRS_struct.p.Navg(ii), MRS_struct.p.npoints(ii));
            PRIAMWaterData = zeros(length(MRS_struct.p.vox), MRS_struct.p.Nwateravg(ii), MRS_struct.p.npoints(ii));
            for kk = 1:MRS_struct.p.Navg(ii)
                PRIAMData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(MRS_struct.fids.data(:,kk,:));
                % Phase by multiplying with normalized complex conjugate of first point
                conj_norm = conj(PRIAMData(:,kk,1)) ./ abs(conj(PRIAMData(:,kk,1)));
                PRIAMData(:,kk,:) = PRIAMData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints(ii)]);
            end
            for kk = 1:MRS_struct.p.Nwateravg(ii)
                PRIAMWaterData(:,kk,:) = MRS_struct.p.SENSE.U * squeeze(MRS_struct.fids.data_water(:,kk,:));
                % Phase by multiplying with normalized complex conjugate of first point
                conj_norm = conj(PRIAMWaterData(:,kk,1)) ./ abs(conj(PRIAMWaterData(:,kk,1)));
                PRIAMWaterData(:,kk,:) = PRIAMWaterData(:,kk,:) .* repmat(conj_norm, [1 1 MRS_struct.p.npoints(ii)]);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   5. Apply appropriate pre-processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for kk = 1:length(vox) % loop over number of voxels
            
            % Select data from first voxel
            if MRS_struct.p.PRIAM
                MRS_struct.fids.data = squeeze(-PRIAMData(kk,:,:))';
                MRS_struct.fids.data_water = squeeze(PRIAMWaterData(kk,:,:))';
            end
            
            % Zero-fill to obtain nominal spectral resolution of 0.061 Hz/point
            MRS_struct.p.ZeroFillTo(ii) = round(32768 / 2000 * MRS_struct.p.sw(ii));
            MRS_struct.p.zf(ii) = MRS_struct.p.ZeroFillTo(ii) / MRS_struct.p.npoints(ii);
            time = (1:size(MRS_struct.fids.data,1)) / MRS_struct.p.sw(ii);
            
            % Finish processing water data
            if strcmp(MRS_struct.p.reference, 'H2O')
                
                % Some instances where the water signal has not been averaged yet
                if MRS_struct.p.join
                    MRS_struct.fids.data_water = mean(MRS_struct.fids.data_water,2);
                elseif strcmp(MRS_struct.p.vendor, 'Philips_raw')
                    MRS_struct.fids.data_water = mean(MRS_struct.fids.data_water(kk,:,:),2);
                end
                
                % Performing eddy current corrrection on the water-suppressed data
                if MRS_struct.p.metab_ECC
                    MRS_struct.fids.data = EddyCurrentCorrection(MRS_struct.fids.data, MRS_struct.fids.data_water);
                end
                
                % Performing eddy current corrrection on the unsuppressed water data
                if MRS_struct.p.water_ECC
                    MRS_struct.fids.data_water = EddyCurrentCorrection(MRS_struct.fids.data_water, MRS_struct.fids.data_water);
                end
                
                % Line-broadening, zero-filling and FFT
                % Water data may have different bandwidth
                if isfield(MRS_struct.p, 'sw_water')
                    time_water = (1:size(MRS_struct.fids.data_water,1)) / MRS_struct.p.sw_water(ii);
                else
                    time_water = (1:size(MRS_struct.fids.data_water,1)) / MRS_struct.p.sw(ii);
                end
                MRS_struct.fids.data_water = MRS_struct.fids.data_water .* exp(-time_water' * MRS_struct.p.LB * pi);
                MRS_struct.spec.(vox{kk}).water(ii,:) = fftshift(fft(MRS_struct.fids.data_water, MRS_struct.p.ZeroFillTo(ii), 1),1);
                
            end
            
            % Global zero-order phase correction
            if ~MRS_struct.p.phantom
                MRS_struct.fids.data = PhaseCorrection(MRS_struct.fids.data, MRS_struct);
            end
            
            % Line-broadening, zero-filling and FFT
            AllFramesFT = MRS_struct.fids.data .* repmat(exp(-time' * MRS_struct.p.LB * pi), [1 size(MRS_struct.fids.data,2)]);
            AllFramesFT = fftshift(fft(AllFramesFT, MRS_struct.p.ZeroFillTo(ii), 1),1);
            
            % Work out ppm axis
            freqRange = MRS_struct.p.sw(ii) / MRS_struct.p.LarmorFreq(ii);
            if MRS_struct.p.phantom
                F0 = 4.8;
            else
                F0 = 4.68;
            end
            MRS_struct.spec.freq = (MRS_struct.p.ZeroFillTo(ii) + 1 - (1:1:MRS_struct.p.ZeroFillTo(ii))) / MRS_struct.p.ZeroFillTo(ii) * freqRange + F0 - freqRange/2;
            
            MRS_struct.p.dt(ii)             = 1/MRS_struct.p.sw(ii);
            MRS_struct.p.SpecRes(ii)        = MRS_struct.p.sw(ii) / MRS_struct.p.npoints(ii);
            MRS_struct.p.SpecResNominal(ii) = MRS_struct.p.sw(ii) / MRS_struct.p.ZeroFillTo(ii);
            MRS_struct.p.Tacq(ii)           = 1/MRS_struct.p.SpecRes(ii);
            
            % Frame-by-frame determination of frequency of residual water or Cr (if HERMES/HERCULES or GSH editing)
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.15 & MRS_struct.spec.freq - 3.02 <= 0.15;
            else
                F0freqRange = MRS_struct.spec.freq - F0 >= -0.2 & MRS_struct.spec.freq - F0 <= 0.2;
            end
            [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
            F0freqRange = MRS_struct.spec.freq(F0freqRange);
            MRS_struct.spec.F0freq{ii} = F0freqRange(FrameMaxPos);
            
            % Estimate average amount of F0 offset
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target,'GSH'))
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 3.02);
            elseif any(strcmp(MRS_struct.p.vendor, {'Siemens_rda', 'Siemens_twix', 'Siemens_dicom'}))
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - 4.7); % Siemens assumes 4.7 ppm as F0
            else
                MRS_struct.out.AvgDeltaF0(ii) = mean(F0freqRange(FrameMaxPos) - F0);
            end
            
            % Use frame-by-frame frequency of Cr for RobustSpectralRegistration
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target, 'GSH'))
                MRS_struct.spec.F0freq2{ii} = MRS_struct.spec.F0freq{ii};
            else
                F0freqRange = MRS_struct.spec.freq - 3.02 >= -0.15 & MRS_struct.spec.freq - 3.02 <= 0.15;
                [~,FrameMaxPos] = max(abs(real(AllFramesFT(F0freqRange,:))),[],1);
                F0freqRange = MRS_struct.spec.freq(F0freqRange);
                MRS_struct.spec.F0freq2{ii} = F0freqRange(FrameMaxPos);
            end
            
            % Frame-by-frame alignment
            switch MRS_struct.p.alignment
                case {'Cr','Cho','NAA'}
                    [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFT, MRS_struct);
                case 'H2O'
                    [AllFramesFTrealign, MRS_struct] = AlignUsingH2O(AllFramesFT, MRS_struct);
                case 'SpecReg'
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistration(MRS_struct,0);
                case 'SpecRegDual'
                    % Dual-channel spectral registration is applied separately to ON and OFF and they are coregistered after
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistration(MRS_struct,0,1);
                case 'SpecRegHERMES'
                    [AllFramesFTrealign, MRS_struct] = SpectralRegistrationHERMES(MRS_struct);
                case 'RobustSpecReg'
                    [AllFramesFTrealign, MRS_struct] = RobustSpectralRegistration(MRS_struct);
                case 'none'
                    % do nothing
                    AllFramesFTrealign = AllFramesFT;
                    MRS_struct.out.reject{ii} = zeros(1,size(AllFramesFT,2));
                otherwise
                    error('Alignment method in GannetPreInitialise.m not recognized. Check spelling.');
            end
            
            MRS_struct.spec.AllFramesFT        = AllFramesFT;
            MRS_struct.spec.AllFramesFTrealign = AllFramesFTrealign;
            
            % Average subspectra and generate DIFF spectra
            MRS_struct = SignalAveraging(MRS_struct, AllFramesFT, AllFramesFTrealign, ii, kk, vox);
            
            % Remove residual water from diff and diff_noalign spectra using HSVD
            if MRS_struct.p.water_removal
                
                for jj = 1:length(MRS_struct.p.target)
                    if jj == 1
                        fprintf('\nRemoving the residual water signal using HSVD...\n');
                    end
                    
                    % Convert DIFF spectra to time domain, apply water filter, convert back to frequency domain
                    fids.diff = WaterRemovalHSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = fftshift(fft(fids.diff));
                    
                    fids.diff_noalign = WaterRemovalHSVD(ifft(ifftshift(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:).')), ...
                        MRS_struct.p.sw(ii)/1e3, 8, -0.08, 0.08, 0, 2048);
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = fftshift(fft(fids.diff_noalign));
                    
                    % Need to perform baseline correction on filtered data
                    freqbounds           = MRS_struct.spec.freq <= 8 & MRS_struct.spec.freq >= 7;
                    baseMean_diff        = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,freqbounds)));
                    baseMean_diffnoalign = mean(real(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,freqbounds)));
                    
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) - baseMean_diff;
                    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) - baseMean_diffnoalign;
                end

            else

                fprintf('\n');
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   6. Build GannetLoad output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ishandle(101)
                clf(101);
            end
            if MRS_struct.p.hide
                h = figure('Visible', 'off');
            else
                h = figure(101);
            end
            % Open figure in center of screen
            scr_sz = get(0,'ScreenSize');
            fig_w = 1000;
            fig_h = 707;
            set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
            set(h,'Color',[1 1 1]);
            figTitle = 'GannetLoad Output';
            set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
            
            % Top left
            if length(MRS_struct.p.target) == 3
                subplot(5,2,1:2:5);
            else
                subplot(2,2,1);
            end
            PlotPrePostAlign(MRS_struct, vox, ii, kk);
            
            % Top right
            if MRS_struct.p.phantom
                if MRS_struct.p.HERMES
                    F0 = 3.02;
                else
                    F0 = 4.8;
                end
            elseif MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target, 'GSH'))
                F0 = 3.02;
            else
                F0 = 4.68;
            end
            
            subplot(2,2,2);
            if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                rejectframesplot = (1./MRS_struct.out.reject{ii}) .* MRS_struct.spec.F0freq{ii};
            end
            hold on;
            plot([1 size(MRS_struct.fids.data,2)], [F0 F0], '-k')
            plot([1 size(MRS_struct.fids.data,2)], [F0-0.04 F0-0.04], '--k')
            plot([1 size(MRS_struct.fids.data,2)], [F0+0.04 F0+0.04], '--k');
            plot(1:size(MRS_struct.fids.data,2), MRS_struct.spec.F0freq{ii}', 'Color', 'b');
            if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                plot(1:size(MRS_struct.fids.data,2), rejectframesplot, 'ro');
            end
            hold off;
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target, 'GSH'))
                text(size(MRS_struct.fids.data,2) + 0.025*size(MRS_struct.fids.data,2), F0, {'Nominal','Cr freq.'}, 'FontSize', 8);
            else
                text(size(MRS_struct.fids.data,2) + 0.025*size(MRS_struct.fids.data,2), F0, {'Nominal','water freq.'}, 'FontSize', 8);
            end
            set(gca,'TickDir','out','box','off','XLim',[1 size(MRS_struct.fids.data,2)], ...
                'YLim',[min([F0-0.06 MRS_struct.spec.F0freq{ii}-0.005]) max([F0+0.06 MRS_struct.spec.F0freq{ii}+0.005])]);
            if size(MRS_struct.fids.data,2) == 2
                set(gca,'XTick',[1 2]);
            end
            xlabel('average');
            ylabel('ppm');
            if MRS_struct.p.HERMES || any(strcmp(MRS_struct.p.target, 'GSH'))
                title('Cr frequency');
            else
                title('Water frequency');
            end
            
            % Bottom left
            if length(MRS_struct.p.target) == 3
                subplot(3,2,5);
            else
                subplot(2,2,3);
            end
            if ~strcmp(MRS_struct.p.alignment, 'none')
                CrFitLimLow = 2.72;
                CrFitLimHigh = 3.12;
                plotrange = MRS_struct.spec.freq <= CrFitLimHigh & MRS_struct.spec.freq >= CrFitLimLow;
                CrFitRange = sum(plotrange);
                plotrealign = [real(AllFramesFT(plotrange,:)); real(AllFramesFTrealign(plotrange,:))];
                % Don't display rejects
                if ~MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4
                    plotrealign(CrFitRange+1:end,(MRS_struct.out.reject{ii} == 1)) = min(plotrealign(:));
                end
                imagesc(plotrealign);
                colormap('parula');
                title({'Cr frequency','(pre- and post-alignment)'});
                xlabel('average');
                ylabel('ppm');
                set(gca, 'YTick', [CrFitRange * (CrFitLimHigh - 3.02) / (CrFitLimHigh - CrFitLimLow) ...
                                   CrFitRange ...
                                   CrFitRange + CrFitRange * (CrFitLimHigh - 3.02) / (CrFitLimHigh - CrFitLimLow) ...
                                   CrFitRange * 2], ...
                    'YTickLabel', [3.02 CrFitLimLow 3.02 CrFitLimLow], ...
                    'XLim', [1 size(MRS_struct.fids.data,2)], ...
                    'YLim', [1 CrFitRange * 2], ...
                    'TickDir','out','box','off');
                if size(MRS_struct.fids.data,2) == 2
                    set(gca, 'XTick', [1 2]);
                end
                % Add in labels for pre/post
                text(size(plotrealign,2)/18*17, 0.4*size(plotrealign,1), 'PRE', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
                text(size(plotrealign,2)/18*17, 0.9*size(plotrealign,1), 'POST', 'Color', [1 1 1], 'HorizontalAlignment', 'right');
            else
                CrFitLimLow = 2.72;
                CrFitLimHigh = 3.12;
                plotrange = MRS_struct.spec.freq <= CrFitLimHigh & MRS_struct.spec.freq >= CrFitLimLow;
                CrFitRange = sum(plotrange);
                plotrealign = real(AllFramesFTrealign(plotrange,:));
                imagesc(plotrealign);
                colormap('parula');
                title({'Cr frequency','(no alignment)'});
                xlabel('average');
                ylabel('ppm');
                set(gca, 'YTick', [CrFitRange * (CrFitLimHigh - 3.02) / (CrFitLimHigh - CrFitLimLow) ...
                                   CrFitRange], ...
                    'YTickLabel', [3.02 CrFitLimLow], ...
                    'XLim', [1 size(MRS_struct.fids.data,2)], ...
                    'YLim', [1 CrFitRange], ...
                    'TickDir','out','box','off');
                if size(MRS_struct.fids.data,2) == 2
                    set(gca, 'XTick', [1 2]);
                end
            end
            
            % Bottom right
            subplot(2,2,4);
            axis off;
            
            if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii*2-1});
            else
                [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii});
            end
            fname = [tmp tmp2];
            if length(fname) > 30
                fname = sprintf([fname(1:floor((end-1)/2)) '...\n     ' fname(ceil(end/2):end)]);
                shift = 0.02;
            else
                shift = 0;
            end
            text(0.25, 1, 'Filename: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            if MRS_struct.p.join
                text(0.275, 1+shift, [fname ' (+ ' num2str(MRS_struct.p.numFilesPerScan - 1) ' more)'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            else
                text(0.275, 1+shift, fname, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            end
            
            vendor = MRS_struct.p.vendor;
            ind = strfind(vendor,'_');
            if ~isempty(ind)
                vendor(ind:end) = '';
            end
            text(0.25, 0.9, 'Vendor: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.9, vendor, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            
            text(0.25, 0.8, 'TE/TR: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.8, [num2str(MRS_struct.p.TE(ii)) '/' num2str(MRS_struct.p.TR(ii)) ' ms'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            
            text(0.25, 0.7, 'Averages: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.7, num2str(MRS_struct.p.Navg(ii)), 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            
            tmp = [num2str(MRS_struct.p.voxdim(ii,1)) ' \times ' num2str(MRS_struct.p.voxdim(ii,2)) ' \times ' num2str(MRS_struct.p.voxdim(ii,3)) ' mm^{3}'];
            text(0.25, 0.6, 'Volume: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.6, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            
            text(0.25, 0.5, 'Spectral width: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.5, [num2str(MRS_struct.p.sw(ii)) ' Hz'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
            
            text(0.25, 0.4, 'Data points: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.4, [num2str(MRS_struct.p.npoints(ii)) ' (zero-filled to ' num2str(MRS_struct.p.ZeroFillTo(ii)) ')'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
                        
            text(0.25, 0.3, 'Alignment: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            if strcmp(MRS_struct.p.alignment, 'RobustSpecReg') && MRS_struct.p.use_prealign_ref
                text(0.275, 0.3, [MRS_struct.p.alignment ' (PreAlignRef)'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            else
                text(0.275, 0.3, MRS_struct.p.alignment, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            end
            
            tmp = [num2str(MRS_struct.p.LB) ' Hz'];
            text(0.25, 0.2, 'Line-broadening: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0.2, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            
            text(0.25, 0.1, 'Rejects: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            if MRS_struct.p.weighted_averaging
                text(0.275, 0.1, 'n/a - wgt. avg. used', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            else
                text(0.275, 0.1, num2str(sum(MRS_struct.out.reject{ii})), 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            end
            
            text(0.25, 0, 'LoadVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'HorizontalAlignment', 'right');
            text(0.275, 0, MRS_struct.version.load, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);
            
            % Save output as PDF
            run_count = SavePDF(h, MRS_struct, ii, 1, kk, vox, mfilename, run_count);
            
            % Reorder structure
            if isfield(MRS_struct, 'waterfile')
                structorder = {'version', 'ii', 'metabfile', ...
                    'waterfile', 'p', 'fids', 'spec', 'out'};
            else
                structorder = {'version', 'ii', 'metabfile', ...
                    'p', 'fids', 'spec', 'out'};
            end
            MRS_struct = orderfields(MRS_struct, structorder);
            
        end % end of output loop over voxels
        
    catch ME
        
        fprintf('\n');
        warning('********** An error occured while loading dataset: ''%s''. Check data. Skipping to next dataset in batch. **********', MRS_struct.metabfile{1,ii});
        error_report{catch_ind} = strrep(sprintf(['Filename: %s\n\n' getReport(ME,'extended','hyperlinks','off') ...
            '\n\nVisit https://markmikkelsen.github.io/Gannet-docs/index.html for help.'], MRS_struct.metabfile{1,ii}), '\', '\\');
        catch_ind = catch_ind + 1;
        
    end % end of load-and-processing loop over datasets
    
    % Display report if errors occurred
    if ~isempty(error_report{1}) && ii == MRS_struct.p.numScans
        opts = struct('WindowStyle', 'non-modal', 'Interpreter', 'tex');
        for ll = flip(1:size(error_report,2))
            errordlg(['\fontsize{13}' regexprep(error_report{ll}, '_', '\\_')], sprintf('GannetLoad Error Report (%d of %d)', ll, size(error_report,2)), opts);
        end
    end
    
end

if MRS_struct.p.mat % save MRS_struct as mat file
    mat_name = fullfile(pwd, ['MRS_struct_' vox{kk} '.mat']);
    if exist(mat_name, 'file')
        fprintf('\nUpdating results in %s\n', ['MRS_struct_' vox{kk} '.mat...']);
    else
        fprintf('\nSaving results to %s\n', ['MRS_struct_' vox{kk} '.mat...']);
    end
    save(mat_name, 'MRS_struct', '-v7.3');
end

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','MATLAB:rankDeficientMatrix');

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end



