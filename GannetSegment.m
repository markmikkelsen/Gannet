function MRS_struct = GannetSegment(MRS_struct)
% Relies on SPM12 being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

if ~isstruct(MRS_struct)
    fprintf('\n');
    error('The first input argument ''%s'' must be a structure.', MRS_struct);
end

MRS_struct.info.datetime.segment = datetime('now');
MRS_struct.info.version.segment = '250911';

warning('off'); % temporarily suppress warning messages

% First check if SPM12 is installed and on the search path
spm_version = fileparts(which('spm'));
if isempty(spm_version)
    msg = 'SPM not found! Please install SPM12 and make sure it is in your search path.';
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12/', 'SPM12', msg);
    fprintf('\n');
    error(msg);
elseif strcmpi(spm_version(end-3:end), 'spm8')
    msg = ['SPM8 detected. Gannet no longer supports SPM8. ' ...
           'Please install SPM12 and make sure it is in your search path.'];
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12/', 'SPM12', msg);
    fprintf('\n');
    error(msg);
end

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end

run_count = 0;
setup_spm = 1;

% Loop over voxels if PRIAM
for kk = 1:length(vox)

    for ii = 1:MRS_struct.p.numScans

        % 1. Take NIfTI from GannetCoRegister and segment it in SPM

        [T1dir, T1name, T1ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        struc = MRS_struct.mask.(vox{kk}).T1image{ii};

        % Check to see if segmentation has already been done (and all
        % probability tissue maps are present)
        if ~MRS_struct.p.bids
            tissue_maps = {[T1dir filesep 'c1' T1name T1ext]
                           [T1dir filesep 'c2' T1name T1ext]
                           [T1dir filesep 'c3' T1name T1ext]
                           [T1dir filesep 'c6' T1name T1ext]};
            filesExist = zeros(length(tissue_maps),1);
            for jj = 1:length(tissue_maps)
                filesExist(jj) = exist(tissue_maps{jj}, 'file');
            end
        else % BIDSify
            % Find BIDS probabilistic tissue maps, if there are any
            bids_file = bids.File(MRS_struct.mask.(vox{kk}).T1image{ii});
            if ~exist(fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path), 'dir')
                bids.util.mkdir(fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path));
            end
            input = mergestructs(bids_file.entities, struct('space', 'orig'));
            bids_file.entities = input;
            bids_file.suffix = 'probseg';
            tiss_class = {'GM','WM','CSF','BG'};
            probseg_fname = cell(length(tiss_class),1);
            filesExist = zeros(length(tiss_class),1);
            for jj = 1:length(tiss_class)
                bids_file.entities = mergestructs(input, struct('label', tiss_class{jj}));
                probseg_fname{jj} = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
                filesExist(jj) = exist(probseg_fname{jj}, 'file');
            end
        end

        files_segmented = 0;
        if ~all(filesExist)
            if setup_spm
                % Set up SPM for batch processing (do it once per batch)
                spm('defaults','fmri');
                spm_jobman('initcfg');
                setup_spm = 0;
            end
            if kk == 1 && ii == 1
                fprintf('\nSegmenting %s...', [T1name T1ext]);
            else
                fprintf('Segmenting %s...', [T1name T1ext]);
            end
            CallSPM12segmentation(struc);
        else
            files_segmented = 1;
            if kk == 1 && ii == 1
                fprintf('\n%s has already been segmented...\n', [T1name T1ext]);
            else
                fprintf('%s has already been segmented...\n', [T1name T1ext]);
            end
        end
        
        % 2. Calculate QC metrics and GM, WM, and CSF fractions for each voxel
        
        if strcmp(T1dir, '')
            T1dir = '.';
        end

        % Tissue ﻿probability maps
        GM  = [T1dir filesep 'c1' T1name T1ext];
        WM  = [T1dir filesep 'c2' T1name T1ext];
        CSF = [T1dir filesep 'c3' T1name T1ext];
        BG  = [T1dir filesep 'c6' T1name T1ext];

        % Forward deformation field
        [struc_dir, struc_name, struc_ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        MRS_struct.mask.(vox{kk}).fwd_def{ii,:} = fullfile(struc_dir, ['y_' struc_name struc_ext]);

        % BIDSify
        if MRS_struct.p.bids
            seg_mat = fullfile(struc_dir, [struc_name '_seg8.mat']);
            if exist(seg_mat, 'file')
                delete(seg_mat);
            end

            if ~files_segmented
                movefile(GM, probseg_fname{1});
                movefile(WM, probseg_fname{2});
                movefile(CSF, probseg_fname{3});
                movefile(BG, probseg_fname{4});
            end

            bids_file = bids.File(MRS_struct.mask.(vox{kk}).T1image{ii});
            input = mergestructs(bids_file.entities, struct('desc', 'fwddef'));
            bids_file.entities = input;

            fwd_def = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
            if ~files_segmented
                movefile(MRS_struct.mask.(vox{kk}).fwd_def{ii,:}, fwd_def);
            end
            MRS_struct.mask.(vox{kk}).fwd_def{ii,:} = fwd_def;

            GM  = probseg_fname{1};
            WM  = probseg_fname{2};
            CSF = probseg_fname{3};
            BG  = probseg_fname{4};
        end

        GM_vol  = spm_vol(GM);
        WM_vol  = spm_vol(WM);
        CSF_vol = spm_vol(CSF);
        BG_vol  = spm_vol(BG);

        % Segmentation quality metrics (Chua et al. JMRI, 2009; Ganzetti et
        % al. Front. Neuroinform., 2016; Esteban et al. PLOS One, 2017)
        T1     = spm_vol(struc);
        T1_tmp = T1.private.dat(:,:,:);
        
        WM_vol_tmp = WM_vol.private.dat(:,:,:);
        WM_vol_tmp(WM_vol_tmp < 0.9) = NaN; % threshold at 0.9 to avoid partial volume effects
        WM_vol_thresh = WM_vol_tmp .* T1_tmp;
        WM_vol_thresh = WM_vol_thresh(:);
        
        GM_vol_tmp = GM_vol.private.dat(:,:,:);
        GM_vol_tmp(GM_vol_tmp < 0.9) = NaN;
        GM_vol_thresh = GM_vol_tmp .* T1_tmp;
        GM_vol_thresh = GM_vol_thresh(:);
        
        BG_vol_tmp = BG_vol.private.dat(:,:,:);
        BG_vol_tmp(BG_vol_tmp < 0.9) = NaN;
        BG_vol_thresh = BG_vol_tmp .* T1_tmp;
        BG_vol_thresh = BG_vol_thresh(:);
        
        MRS_struct.out.QA.CV_WM(ii) = std(WM_vol_thresh, 'omitnan') / mean(WM_vol_thresh, 'omitnan');
        MRS_struct.out.QA.CV_GM(ii) = std(GM_vol_thresh, 'omitnan') / mean(GM_vol_thresh, 'omitnan');
        MRS_struct.out.QA.CJV(ii)   = (std(WM_vol_thresh, 'omitnan') + std(GM_vol_thresh, 'omitnan')) ...
                                          / abs(mean(WM_vol_thresh, 'omitnan') - mean(GM_vol_thresh, 'omitnan'));
        MRS_struct.out.QA.SNR(ii)   = std([WM_vol_thresh; GM_vol_thresh], 'omitnan') / mean([WM_vol_thresh; GM_vol_thresh], 'omitnan');
        MRS_struct.out.QA.CNR(ii)   = abs(mean(WM_vol_thresh, 'omitnan') - mean(GM_vol_thresh, 'omitnan')) / ...
                                          sqrt(var(BG_vol_thresh, 'omitnan') + var(WM_vol_thresh, 'omitnan') + var(GM_vol_thresh, 'omitnan'));
        
        T1_tmp  = T1_tmp(:);
        n_vox   = numel(T1_tmp);
        efc_max = n_vox * (1/sqrt(n_vox)) * log(1/sqrt(n_vox));
        b_max   = sqrt(sum(T1_tmp.^2));
        MRS_struct.out.QA.EFC(ii) = (1/efc_max) .* sum((T1_tmp ./ b_max) .* log((T1_tmp + eps) ./ b_max));
        
        % Voxel mask
        vox_mask_vol = spm_vol(MRS_struct.mask.(vox{kk}).fname{ii});
        [a,b,c] = fileparts(vox_mask_vol.fname);
        
        % GM
        if MRS_struct.p.bids
            bids_file = bids.File(MRS_struct.mask.(vox{kk}).fname{ii});
            input = mergestructs(bids_file.entities, struct('label', 'GM'));
            bids_file.entities = input;
            GM_vox.fname = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
        else
            GM_vox.fname = fullfile(a, [b '_GM' c]);
        end
        GM_vox.descrip = 'MRS_voxel_mask_GM';
        GM_vox.dim = vox_mask_vol.dim;
        GM_vox.dt = vox_mask_vol.dt;
        GM_vox.mat = vox_mask_vol.mat;
        GM_vox_mask_vol = GM_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        GM_vox = spm_write_vol(GM_vox, GM_vox_mask_vol);
        
        % WM
        if MRS_struct.p.bids
            bids_file = bids.File(MRS_struct.mask.(vox{kk}).fname{ii});
            input = mergestructs(bids_file.entities, struct('label', 'WM'));
            bids_file.entities = input;
            WM_vox.fname = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
        else
            WM_vox.fname = fullfile(a, [b '_WM' c]);
        end
        WM_vox.descrip = 'MRS_voxel_mask_WM';
        WM_vox.dim = vox_mask_vol.dim;
        WM_vox.dt = vox_mask_vol.dt;
        WM_vox.mat = vox_mask_vol.mat;
        WM_voxmask_vol = WM_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        WM_vox = spm_write_vol(WM_vox, WM_voxmask_vol);
        
        % CSF
        if MRS_struct.p.bids
            bids_file = bids.File(MRS_struct.mask.(vox{kk}).fname{ii});
            input = mergestructs(bids_file.entities, struct('label', 'CSF'));
            bids_file.entities = input;
            CSF_vox.fname = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
        else
            CSF_vox.fname = fullfile(a, [b '_CSF' c]);
        end
        CSF_vox.descrip = 'MRS_voxel_mask_CSF';
        CSF_vox.dim = vox_mask_vol.dim;
        CSF_vox.dt = vox_mask_vol.dt;
        CSF_vox.mat = vox_mask_vol.mat;
        CSF_voxmask_vol = CSF_vol.private.dat(:,:,:) .* vox_mask_vol.private.dat(:,:,:);
        CSF_vox = spm_write_vol(CSF_vox, CSF_voxmask_vol);
        
        % 3. Calculate a CSF-corrected i.u. value and output it to the structure
        
        GM_vox_n  = GM_vox.private.dat(:,:,:);
        GM_sum    = sum(GM_vox_n(GM_vox_n > 0.9)); % threshold at 0.9 to improve accuracy
        WM_vox_n  = WM_vox.private.dat(:,:,:);
        WM_sum    = sum(WM_vox_n(WM_vox_n > 0.9));
        CSF_vox_n = CSF_vox.private.dat(:,:,:);
        CSF_sum   = sum(CSF_vox_n(CSF_vox_n > 0.9));
        
        fGM  = GM_sum / (GM_sum + WM_sum + CSF_sum);
        fWM  = WM_sum / (GM_sum + WM_sum + CSF_sum);
        fCSF = CSF_sum / (GM_sum + WM_sum + CSF_sum);
        
        MRS_struct.out.(vox{kk}).tissue.fGM(ii)  = fGM;
        MRS_struct.out.(vox{kk}).tissue.fWM(ii)  = fWM;
        MRS_struct.out.(vox{kk}).tissue.fCSF(ii) = fCSF;
        
        % Correction of institutional units only feasible if water-scaling
        % is performed, skip otherwise
        if strcmp(MRS_struct.p.reference, 'H2O')
            target = [MRS_struct.p.target, {'Cr'}, {'Cho'}, {'NAA'}, {'Glu'}]; % Add Cr, Cho, NAA, and Glu
            for jj = 1:length(target)
                if strcmp(target{jj}, 'GABAGlx')
                    MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).GABA.ConcIU(ii) / (1 - fCSF);
                    MRS_struct.out.(vox{kk}).Glx.ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).Glx.ConcIU(ii) / (1 - fCSF);
                else
                    MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_CSFcorr(ii) = ...
                        MRS_struct.out.(vox{kk}).(target{jj}).ConcIU(ii) / (1 - fCSF);
                end
            end
        end

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
        
        if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
            [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii});
        end
        fname = [name ext];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, text_pos-0.12, 'Filename: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, [' ' fname], 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        [~,name,ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [name ext];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        text(0.5, text_pos-0.24, 'Anatomical image: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, [' ' T1image], 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        % Print correction of institutional units only feasible if water-scaling is performed, skip otherwise
        text_pos = text_pos-0.24;
        if strcmp(MRS_struct.p.reference, 'H2O')
            target = MRS_struct.p.target;
            is_GABAGlx = strcmp(target, 'GABAGlx');
            if any(is_GABAGlx)
                if MRS_struct.p.HERMES
                    target = {'GABA','Glx',target{~is_GABAGlx}};
                else
                    target = {'GABA','Glx'};
                end
            end
            for jj = 1:length(target)
                text_pos = text_pos - 0.12;
                switch target{jj}
                    case 'GABA'
                        str1 = 'GABA+/Water (CSF-corrected): ';
                        str2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii));
                    case 'Lac'
                        str1 = 'Lac+MM/Water (CSF-corrected): ';
                        str2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).Lac.ConcIU_CSFcorr(ii));
                    case {'Glx','GSH','EtOH'}
                        str1 = [target{jj} '/Water (CSF-corrected): '];
                        str2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_CSFcorr(ii));
                end
                text(0.5, text_pos, str1, 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
                text(0.5, text_pos, str2, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
            end
        end
        
        str = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fGM(ii));
        text(0.5, text_pos-0.12, 'GM voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, str, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        str = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fWM(ii));
        text(0.5, text_pos-0.24, 'WM voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, str, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        str = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fCSF(ii));
        text(0.5, text_pos-0.36, 'CSF voxel fraction: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.36, str, 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        text(0.5, text_pos-0.48, 'SegmentVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.48, [' ' MRS_struct.info.version.segment], 'Units', 'normalized', 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        % Voxel segmentation
        if isfield(MRS_struct.p,'TablePosition')
            voxoff = MRS_struct.p.voxoff(ii,:) + MRS_struct.p.TablePosition(ii,:);
        else
            voxoff = MRS_struct.p.voxoff(ii,:);
        end
        
        % Plot segmented voxel as a montage
        MRS_struct.mask.(vox{kk}).img_montage{ii} = PlotSegmentedVoxels(struc, voxoff, vox_mask_vol, GM_vox, WM_vox, CSF_vox);

        if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
            [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii});
        end
        fname = [name ext];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        [~,name,ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [name ext];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        t = ['Voxel from ' fname ' on ' T1image];
        title(t, 'FontName', 'Arial', 'FontSize', 15, 'Interpreter', 'none');
        
        % Save output as PDF
        run_count = SavePDF(h, MRS_struct, ii, 1, kk, vox, mfilename, run_count);
        
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

    if MRS_struct.p.csv % export MRS_struct fields into csv file
        MRS_struct = ExportToCSV(MRS_struct, vox{kk}, 'segment');
    end

end

warning('on'); % turn warnings back on

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end


function img_montage = PlotSegmentedVoxels(struc, voxoff, voxmaskvol, GM_vox, WM_vox, CSF_vox)

img_t      = flipud(voxel2world_space(spm_vol(struc), voxoff));
mask_t     = flipud(voxel2world_space(voxmaskvol, voxoff));
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



