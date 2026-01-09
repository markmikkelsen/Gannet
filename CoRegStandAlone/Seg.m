function MRS_struct = Seg(MRS_struct)

% Relies on SPM12 being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation
%
% This script does not require any information from GannetFit.
% This is useful if only the tissue segmentation information is supposed to
% be obtained.

loadFile = which('GannetSegment');
fileID = fopen(loadFile, 'rt');
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.info.version.segment = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
MRS_struct.info.version.segment = out.version;

warning('off'); % temporarily suppress warning messages

% First check if SPM12 is installed and on the search path
spm_version = fileparts(which('spm'));
if isempty(spm_version)
    msg = 'SPM not found! Please install SPM12 and make sure it is in your search path.';
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12/', 'SPM12', msg);
    fprintf('\n');
    error(msg);
elseif strcmpi(spm_version(end-3:end),'spm8')
    msg = ['SPM8 detected. Gannet no longer supports SPM8. ' ...
           'Please install SPM12 and make sure it is in your search path.'];
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12/', 'SPM12', msg);
    fprintf('\n');
    error(msg);
end

vox = MRS_struct.p.vox(1);
kk = 1;
setup_spm = 1;
prob_threshold = 0.9; % threshold to reduce partial volume effect and
                      % improve accuracy of tissue probability maps

for ii = 1:length(MRS_struct.metabfile)
    
    % 1. Take NIfTI from GannetCoRegister and segment it in SPM
    
    [T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
    struc = MRS_struct.mask.(vox{kk}).T1image{ii};
    
    % Check to see if segmentation has already been done (and all
    % probability tissue maps are present)
    tissue_maps = {[T1dir filesep 'c1' T1name T1ext]
                   [T1dir filesep 'c2' T1name T1ext]
                   [T1dir filesep 'c3' T1name T1ext]
                   [T1dir filesep 'c6' T1name T1ext]};
    filesExist = zeros(1,length(tissue_maps));
    for jj = 1:length(tissue_maps)
        filesExist(jj) = exist(tissue_maps{jj}, 'file');
    end
    if ~all(filesExist)
        if setup_spm
            % Set up SPM for batch processing (do it once per batch)
            spm('defaults','fmri');
            spm_jobman('initcfg');
            setup_spm = 0;
        end
        if ii == 1
            fprintf('\nSegmenting %s...', [T1name T1ext]);
        else
            fprintf('Segmenting %s...', [T1name T1ext]);
        end
        CallSPM12segmentation(struc);
    else
        if ii == 1
            fprintf('\n%s has already been segmented...\n', [T1name T1ext]);
        else
            fprintf('%s has already been segmented...\n', [T1name T1ext]);
        end
    end
    
    % 2. Calculate QC metrics and GM, WM, and CSF fractions for each voxel
    
    if strcmp(T1dir,'')
        T1dir = '.';
    end
    
    GM  = [T1dir filesep 'c1' T1name T1ext];
    WM  = [T1dir filesep 'c2' T1name T1ext];
    CSF = [T1dir filesep 'c3' T1name T1ext];
    BG  = [T1dir filesep 'c6' T1name T1ext];

    % Forward deformation field
    [struc_dir, struc_name, struc_ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
    MRS_struct.mask.(vox{kk}).fwd_def{ii,:} = fullfile(struc_dir, ['y_' struc_name struc_ext]);

    GM_vol  = spm_vol(GM);
    WM_vol  = spm_vol(WM);
    CSF_vol = spm_vol(CSF);
    BG_vol  = spm_vol(BG);

    % MRIQC image quality metrics (Esteban et al., 2017,doi:﻿10.1371/journal.pone.0184661;
    % also see: Chua et al., 2009, doi:﻿10.1002/jmri.21768; Ganzetti et al., 2016,
    % doi:﻿10.3389/fninf.2016.00010)
    T1     = spm_vol(struc);
    T1_tmp = T1.private.dat(:,:,:);

    WM_vol_tmp = WM_vol.private.dat(:,:,:);
    WM_vol_tmp(WM_vol_tmp < prob_threshold) = 0;
    T1_WM = WM_vol_tmp .* T1_tmp;
    T1_WM = T1_WM(:);
    T1_WM = T1_WM(T1_WM > 0); % include only nonzero voxels

    GM_vol_tmp = GM_vol.private.dat(:,:,:);
    GM_vol_tmp(GM_vol_tmp < prob_threshold) = 0;
    T1_GM = GM_vol_tmp .* T1_tmp;
    T1_GM = T1_GM(:);
    T1_GM = T1_GM(T1_GM > 0);

    CSF_vol_tmp = CSF_vol.private.dat(:,:,:);
    CSF_vol_tmp(CSF_vol_tmp < prob_threshold) = 0;
    T1_CSF = CSF_vol_tmp .* T1_tmp;
    T1_CSF = T1_CSF(:);
    T1_CSF = T1_CSF(T1_CSF > 0);

    BG_vol_tmp = BG_vol.private.dat(:,:,:);
    BG_vol_tmp(BG_vol_tmp < prob_threshold) = 0;
    T1_BG = BG_vol_tmp .* T1_tmp;
    T1_BG = T1_BG(:);
    T1_BG = T1_BG(T1_BG > 0);

    head_vol_tmp = 1 - BG_vol.private.dat(:,:,:);
    head_vol_tmp(head_vol_tmp < prob_threshold) = 0;
    T1_head = head_vol_tmp .* T1_tmp;
    T1_head = T1_head(:);
    T1_head = T1_head(T1_head > 0);

    MRS_struct.out.QA.CV.WM(ii)  = mad(T1_WM,1) / median(T1_WM);
    MRS_struct.out.QA.CV.GM(ii)  = mad(T1_GM,1) / median(T1_GM);
    MRS_struct.out.QA.CV.CSF(ii) = mad(T1_CSF,1) / median(T1_CSF);
    MRS_struct.out.QA.CJV(ii)    = (mad(T1_WM,1) + mad(T1_GM,1)) / abs(median(T1_WM) - median(T1_GM));
    MRS_struct.out.QA.CNR(ii)    = abs(median(T1_WM) - median(T1_GM)) / sqrt(std(T1_WM).^2 + std(T1_GM).^2 + std(T1_BG).^2);

    T1_tmp  = T1_tmp(:);
    n_vox   = numel(T1_tmp);
    efc_max = n_vox * (1/sqrt(n_vox)) * log(1/sqrt(n_vox));
    b_max   = sqrt(sum(T1_tmp.^2));
    MRS_struct.out.QA.EFC(ii) = (1/efc_max) .* sum((T1_tmp / b_max) .* log((T1_tmp + eps) / b_max));

    MRS_struct.out.QA.FBER(ii)   = median(abs(T1_head).^2) / median(abs(T1_BG).^2);
    MRS_struct.out.QA.WM2MAX(ii) = median(T1_WM) / prctile(T1_tmp, 99.95);

    MRS_struct.out.QA.SNR.WM(ii)    = median(T1_WM) / (std(T1_WM) * sqrt(numel(T1_WM) / (numel(T1_WM) - 1)));
    MRS_struct.out.QA.SNR.GM(ii)    = median(T1_GM) / (std(T1_GM) * sqrt(numel(T1_GM) / (numel(T1_GM) - 1)));
    MRS_struct.out.QA.SNR.CSF(ii)   = median(T1_CSF) / (std(T1_CSF) * sqrt(numel(T1_CSF) / (numel(T1_CSF) - 1)));
    MRS_struct.out.QA.SNR.total(ii) = mean([MRS_struct.out.QA.SNR.WM(ii) MRS_struct.out.QA.SNR.GM(ii) MRS_struct.out.QA.SNR.CSF(ii)]);

    MRS_struct.out.QA.SNR_D.WM(ii)    = median(T1_WM) / (sqrt(2 / (4 - pi)) * mad(T1_BG,1));
    MRS_struct.out.QA.SNR_D.GM(ii)    = median(T1_GM) / (sqrt(2 / (4 - pi)) * mad(T1_BG,1));
    MRS_struct.out.QA.SNR_D.CSF(ii)   = median(T1_CSF) / (sqrt(2 / (4 - pi)) * mad(T1_BG,1));
    MRS_struct.out.QA.SNR_D.total(ii) = mean([MRS_struct.out.QA.SNR_D.WM(ii) MRS_struct.out.QA.SNR_D.GM(ii) MRS_struct.out.QA.SNR_D.CSF(ii)]);

    % Loop over voxels if PRIAM
    for kk = 1:length(vox)

        vox_mask_vol = spm_vol(cell2mat(MRS_struct.mask.(vox{kk}).fname(ii)));
        [a,b,c] = fileparts(vox_mask_vol.fname);

        % GM
        GM_vox.fname = fullfile(a, [b '_GM' c]);
        GM_vox.descrip = 'MRS_voxel_mask_GM';
        GM_vox.dim = vox_mask_vol.dim;
        GM_vox.dt = vox_mask_vol.dt;
        GM_vox.mat = vox_mask_vol.mat;
        GM_vox_mask_vol = GM_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        GM_vox = spm_write_vol(GM_vox, GM_vox_mask_vol);
        
        % WM
        WM_vox.fname = fullfile(a, [b '_WM' c]);
        WM_vox.descrip = 'MRS_voxel_mask_WM';
        WM_vox.dim = vox_mask_vol.dim;
        WM_vox.dt = vox_mask_vol.dt;
        WM_vox.mat = vox_mask_vol.mat;
        WM_vox_mask_vol = WM_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        WM_vox = spm_write_vol(WM_vox, WM_vox_mask_vol);
        
        % CSF
        CSF_vox.fname = fullfile(a, [b '_CSF' c]);
        CSF_vox.descrip = 'MRS_voxel_mask_CSF';
        CSF_vox.dim = vox_mask_vol.dim;
        CSF_vox.dt = vox_mask_vol.dt;
        CSF_vox.mat = vox_mask_vol.mat;
        CSF_vox_mask_vol = CSF_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        CSF_vox = spm_write_vol(CSF_vox, CSF_vox_mask_vol);
        
        % 3. Calculate a CSF-corrected i.u. value and output it to the structure
        
        GM_vox_n  = GM_vox.private.dat(:,:,:);
        GM_sum    = sum(GM_vox_n(GM_vox_n > prob_threshold));
        WM_vox_n  = WM_vox.private.dat(:,:,:);
        WM_sum    = sum(WM_vox_n(WM_vox_n > prob_threshold));
        CSF_vox_n = CSF_vox.private.dat(:,:,:);
        CSF_sum   = sum(CSF_vox_n(CSF_vox_n > prob_threshold));

        fGM  = GM_sum / (GM_sum + WM_sum + CSF_sum);
        fWM  = WM_sum / (GM_sum + WM_sum + CSF_sum);
        fCSF = CSF_sum / (GM_sum + WM_sum + CSF_sum);

        MRS_struct.out.(vox{kk}).tissue.fGM(ii)  = fGM;
        MRS_struct.out.(vox{kk}).tissue.fWM(ii)  = fWM;
        MRS_struct.out.(vox{kk}).tissue.fCSF(ii) = fCSF;

        if MRS_struct.p.normalize
            if setup_spm
                % Set up SPM for batch processing (do it once per batch)
                spm('defaults','fmri');
                spm_jobman('initcfg');
                setup_spm = 0;
            end
            MRS_struct = NormalizeVoxelMask(MRS_struct, vox, ii, kk);
            if kk == length(vox) && ii == MRS_struct.p.numScans && MRS_struct.p.numScans > 1
                MRS_struct = VoxelMaskOverlap(MRS_struct);
            end
        end

        % 4. Build output

        if ishandle(104)
            clf(104);
        end
        if MRS_struct.p.hide
            h = figure('Visible', 'off');
        else
            h = figure(104);
        end
        if ~isMATLABReleaseOlderThan("R2025a")
            h.Theme = 'light';
        end
        % Open figure in center of screen
        scr_sz = get(0,'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetSegment Output';
        set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
        
        % Output results
        ha  = subplot(2,3,4:6);
        pos = get(ha, 'Position');
        set(ha, 'Position', [0 pos(2) 1 pos(4)]);
        axis off;
        
        text_pos = 1;
        
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{ii});
        end
        fname = [tmp2 tmp3];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, text_pos-0.12, 'Filename: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, [' ' fname], 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        [~,tmp1,tmp2] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp1 tmp2];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        text(0.5, text_pos-0.24, 'Anatomical image: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, [' ' T1image], 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fGM(ii));
        text(0.5, text_pos-0.36, 'GM voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.36, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fWM(ii));
        text(0.5, text_pos-0.48, 'WM voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.48, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fCSF(ii));
        text(0.5, text_pos-0.6, 'CSF voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.6, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        text(0.5, text_pos-0.72, 'SegmentVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.72, [' ' MRS_struct.info.version.segment],  'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        if isfield(MRS_struct.p,'TablePosition')
            voxoff = MRS_struct.p.voxoff(ii,:) + MRS_struct.p.TablePosition(ii,:);
        else
            voxoff = MRS_struct.p.voxoff(ii,:);
        end
        
        % Plot segmented voxel as a montage
        MRS_struct.mask.(vox{kk}).img_montage{ii} = PlotSegmentedVoxels(struc, voxoff, vox_mask_vol, GM_vox, WM_vox, CSF_vox);
        
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{ii});
        end
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        [~,tmp3,tmp4] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp3 tmp4];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        t = ['Voxel from ' fname ' on ' T1image];
        title(t, 'FontName', 'Arial', 'FontSize', 15, 'Interpreter', 'none');

        % Gannet logo
        axes('Position', [0.8825, 0.04, 0.125, 0.125], 'Units', 'normalized');
        Gannet_logo = fullfile(fileparts(which('GannetLoad')), 'Gannet3_logo.png');
        I = imread(Gannet_logo, 'BackgroundColor', 'none');
        imshow(I);
        axis off image;

        % Gannet version
        d.left   = 0;
        d.bottom = 0.02;
        d.width  = 1;
        d.height = 0.02;
        axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
        text(0.9925, 0, MRS_struct.info.version.Gannet, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
        axis off;

        % Gannet documentation
        axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
        str = 'For complete documentation, please visit: https://markmikkelsen.github.io/Gannet-docs';
        text(0.5, 0, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11, 'HorizontalAlignment', 'center');
        axis off square;

        % Batch number and output time
        d.bottom = 0.98;
        axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
        text(0.0075, 0, ['Batch file: ' num2str(ii) ' of ' num2str(MRS_struct.p.numScans)], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11, 'HorizontalAlignment', 'left');
        text(0.9925, 0, char(datetime('now','Format','dd-MMM-y HH:mm:ss')), 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11, 'HorizontalAlignment', 'right');
        axis off;

        % For Philips .data
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            fullpath = MRS_struct.metabfile{ii};
            fullpath = regexprep(fullpath, '.data', '_data');
            fullpath = regexprep(fullpath, '\', '_');
            fullpath = regexprep(fullpath, '/', '_');
        end
        
        if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
            [~, metabfile_nopath] = fileparts(MRS_struct.metabfile{ii*2-1});
        else
            [~, metabfile_nopath, ext] = fileparts(MRS_struct.metabfile{ii});
            if strcmpi(ext, '.gz')
                metabfile_nopath(end-3:end) = [];
            end
        end
        
        if any(strcmp(listfonts,'Arial'))
            set(findall(h,'-property','FontName'),'FontName','Arial');
        end
        set(findall(h,'-property','XColor'),'XColor',[0 0 0]);
        set(findall(h,'-property','YColor'),'YColor',[0 0 0]);
        
        % Create output folder
        if ~exist(fullfile(pwd, 'CoRegStandAlone_output'),'dir')
            mkdir(fullfile(pwd, 'CoRegStandAlone_output'));
        end
        
        % Save PDF output
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[11 8.5]);
        set(h,'PaperPosition',[0 0 11 8.5]);
        
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            pdfname = fullfile(pwd, 'CoRegStandAlone_output', [fullpath '_' vox{kk} '_segment.pdf']);
        else
            pdfname = fullfile(pwd, 'CoRegStandAlone_output', [metabfile_nopath '_' vox{kk} '_segment.pdf']);
        end
        saveas(h, pdfname);
        
        
    end
    
end

warning('on'); % turn warnings back on

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide
    close(figTitle);
end


function img_montage = PlotSegmentedVoxels(struc, voxoff, vox_mask_vol, GM_vox, WM_vox, CSF_vox)

img_t      = flipud(voxel2world_space(spm_vol(struc), voxoff));
mask_t     = flipud(voxel2world_space(vox_mask_vol, voxoff));
mask_t_GM  = flipud(voxel2world_space(GM_vox, voxoff));
mask_t_WM  = flipud(voxel2world_space(WM_vox, voxoff));
mask_t_CSF = flipud(voxel2world_space(CSF_vox, voxoff));

w_t = zeros(size(img_t));

tmp       = img_t(:);
img_t     = repmat(img_t / (mean(tmp(tmp > 0.01)) + 3*std(tmp(tmp > 0.01))), [1 1 3]);
img_t_GM  = img_t;
img_t_WM  = img_t;
img_t_CSF = img_t;

c_img_t = zeros(size(img_t));

% Voxel mask
vox_mx = 1;
vox_mn = 0;

mask_t(mask_t(:) < vox_mn) = vox_mn;
mask_t(mask_t(:) > vox_mx) = vox_mx;
mask_t = (mask_t - vox_mn) / (vox_mx - vox_mn);

mask_t_GM(mask_t_GM(:) < vox_mn) = vox_mn;
mask_t_GM(mask_t_GM(:) > vox_mx) = vox_mx;
mask_t_GM = (mask_t_GM - vox_mn) / (vox_mx - vox_mn);

mask_t_WM(mask_t_WM(:) < vox_mn) = vox_mn;
mask_t_WM(mask_t_WM(:) > vox_mx) = vox_mx;
mask_t_WM = (mask_t_WM - vox_mn) / (vox_mx - vox_mn);

mask_t_CSF(mask_t_CSF(:) < vox_mn) = vox_mn;
mask_t_CSF(mask_t_CSF(:) > vox_mx) = vox_mx;
mask_t_CSF = (mask_t_CSF - vox_mn) / (vox_mx - vox_mn);

mask_t     = 0.5 * mask_t;
mask_t_GM  = 0.8 * mask_t_GM;
mask_t_WM  = 0.8 * mask_t_WM;
mask_t_CSF = 0.8 * mask_t_CSF;

vox_color = [1 1 0];

c_mask_t     = c_img_t + cat(3, mask_t * vox_color(1,1), mask_t * vox_color(1,2), mask_t * vox_color(1,3));
c_mask_t_GM  = c_img_t + cat(3, mask_t_GM * vox_color(1,1), mask_t_GM * vox_color(1,2), mask_t_GM * vox_color(1,3));
c_mask_t_WM  = c_img_t + cat(3, mask_t_WM * vox_color(1,1), mask_t_WM * vox_color(1,2), mask_t_WM * vox_color(1,3));
c_mask_t_CSF = c_img_t + cat(3, mask_t_CSF * vox_color(1,1), mask_t_CSF * vox_color(1,2), mask_t_CSF * vox_color(1,3));

w_mask_t     = w_t + mask_t;
w_mask_t_GM  = w_t + mask_t_GM;
w_mask_t_WM  = w_t + mask_t_WM;
w_mask_t_CSF = w_t + mask_t_CSF;

img_t     = repmat(1 - w_mask_t, [1 1 3]) .* img_t + c_mask_t;
img_t_GM  = repmat(1 - w_mask_t_GM, [1 1 3]) .* img_t_GM + c_mask_t_GM;
img_t_WM  = repmat(1 - w_mask_t_WM, [1 1 3]) .* img_t_WM + c_mask_t_WM;
img_t_CSF = repmat(1 - w_mask_t_CSF, [1 1 3]) .* img_t_CSF + c_mask_t_CSF;

img_t(img_t < 0)         = 0; img_t(img_t > 1)         = 1;
img_t_GM(img_t_GM < 0)   = 0; img_t_GM(img_t_GM > 1)   = 1;
img_t_WM(img_t_WM < 0)   = 0; img_t_WM(img_t_WM > 1)   = 1;
img_t_CSF(img_t_CSF < 0) = 0; img_t_CSF(img_t_CSF > 1) = 1;

img_montage = cat(2, img_t, img_t_GM, img_t_WM, img_t_CSF);

hb = subplot(2,3,1:3);
imagesc(img_montage);
axis equal tight off;
text(floor(size(mask_t,2)/2), 20, 'Voxel', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'GM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(2*floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'WM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(3*floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'CSF', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
set(hb, 'Position', [0 0.17 1 1]);



