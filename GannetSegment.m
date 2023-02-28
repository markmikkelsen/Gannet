function MRS_struct = GannetSegment(MRS_struct)

% Relies on SPM12 being installed
%
% Runs segmentation script if segmented images not present according to
% file convention of c1, c2 and c3 as prefixes on the anatomical image name
% for the GM, WM and CSF segmentations. If these files are present, they
% are loaded and used for the voxel segmentation

if nargin == 0
    error('MATLAB:minrhs','Not enough input arguments.');
end

MRS_struct.version.segment = '230228';

warning('off'); % temporarily suppress warning messages

% First check if SPM12 is installed and on the search path
spm_version = fileparts(which('spm'));
if isempty(spm_version)
    msg = 'SPM not found! Please install SPM12 and make sure it is in your search path.';
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12', 'SPM12', msg);
    error(msg);
elseif strcmpi(spm_version(end-3:end),'spm8')
    msg = ['SPM8 detected. Gannet no longer supports SPM8. ' ...
           'Please install SPM12 and make sure it is in your search path.'];
    msg = hyperlink('https://www.fil.ion.ucl.ac.uk/spm/software/spm12', 'SPM12', msg);
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
        tmp = {[T1dir '/c1' T1name T1ext]
               [T1dir '/c2' T1name T1ext]
               [T1dir '/c3' T1name T1ext]
               [T1dir '/c6' T1name T1ext]};
        filesExist = zeros(1,length(tmp));
        for jj = 1:length(tmp)
            filesExist(jj) = exist(tmp{jj}, 'file');
        end
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
            if kk == 1 && ii == 1
                fprintf('\n%s has already been segmented...\n', [T1name T1ext]);
            else
                fprintf('%s has already been segmented...\n', [T1name T1ext]);
            end
        end
        
        % 2. Calculate QC metrics and GM, WM, and CSF fractions for each voxel
        
        if strcmp(T1dir,'')
            T1dir = '.';
        end
        
        % Tissue﻿probability maps
        GM  = [T1dir '/c1' T1name T1ext];
        WM  = [T1dir '/c2' T1name T1ext];
        CSF = [T1dir '/c3' T1name T1ext];
        air = [T1dir '/c6' T1name T1ext];
        
        GMvol  = spm_vol(GM);
        WMvol  = spm_vol(WM);
        CSFvol = spm_vol(CSF);
        airvol = spm_vol(air);
        
        % Segmentation quality metrics (Chua et al. JMRI, 2009; Ganzetti et
        % al. Front. Neuroinform., 2016; Esteban et al. PLOS One, 2017)
        T1     = spm_vol(struc);
        T1_tmp = T1.private.dat(:,:,:);
        
        WMvol_tmp = WMvol.private.dat(:,:,:);
        WMvol_tmp(WMvol_tmp < 0.9) = NaN;
        WMvol_thresh = WMvol_tmp .* T1_tmp;
        WMvol_thresh = WMvol_thresh(:);
        
        GMvol_tmp = GMvol.private.dat(:,:,:);
        GMvol_tmp(GMvol_tmp < 0.9) = NaN;
        GMvol_thresh = GMvol_tmp .* T1_tmp;
        GMvol_thresh = GMvol_thresh(:);
        
        airvol_tmp = airvol.private.dat(:,:,:);
        airvol_tmp(airvol_tmp < 0.9) = NaN;
        airvol_thresh = airvol_tmp .* T1_tmp;
        airvol_thresh = airvol_thresh(:);
        
        MRS_struct.out.QA.CV_WM(ii) = std(WMvol_thresh, 'omitnan') / mean(WMvol_thresh, 'omitnan');
        MRS_struct.out.QA.CV_GM(ii) = std(GMvol_thresh, 'omitnan') / mean(GMvol_thresh, 'omitnan');
        MRS_struct.out.QA.CJV(ii)   = (std(WMvol_thresh, 'omitnan') + std(GMvol_thresh, 'omitnan')) ...
                                          / abs(mean(WMvol_thresh, 'omitnan') - mean(GMvol_thresh, 'omitnan'));
        MRS_struct.out.QA.CNR(ii)   = abs(mean(WMvol_thresh, 'omitnan') - mean(GMvol_thresh, 'omitnan')) / ...
                                          sqrt(var(airvol_thresh, 'omitnan') + var(WMvol_thresh, 'omitnan') + var(GMvol_thresh, 'omitnan'));
        
        T1_tmp  = T1_tmp(:);
        n_vox   = numel(T1_tmp);
        efc_max = n_vox * (1/sqrt(n_vox)) * log(1/sqrt(n_vox));
        b_max   = sqrt(sum(T1_tmp.^2));
        MRS_struct.out.QA.EFC(ii) = (1/efc_max) .* sum((T1_tmp ./ b_max) .* log((T1_tmp + eps) ./ b_max));
        
        % Voxel mask
        voxmaskvol = spm_vol(MRS_struct.mask.(vox{kk}).outfile{ii});
        [a,b,c] = fileparts(voxmaskvol.fname);
        
        % GM
        O_GMvox.fname = fullfile(a, [b '_GM' c]);
        O_GMvox.descrip = 'MRS_voxel_mask_GM';
        O_GMvox.dim = voxmaskvol.dim;
        O_GMvox.dt = voxmaskvol.dt;
        O_GMvox.mat = voxmaskvol.mat;
        GM_voxmask_vol = GMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_GMvox = spm_write_vol(O_GMvox, GM_voxmask_vol);
        
        % WM
        O_WMvox.fname = fullfile(a, [b '_WM' c]);
        O_WMvox.descrip = 'MRS_voxel_mask_WM';
        O_WMvox.dim = voxmaskvol.dim;
        O_WMvox.dt = voxmaskvol.dt;
        O_WMvox.mat = voxmaskvol.mat;
        WM_voxmask_vol = WMvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_WMvox = spm_write_vol(O_WMvox, WM_voxmask_vol);
        
        % CSF
        O_CSFvox.fname = fullfile(a, [b '_CSF' c]);
        O_CSFvox.descrip = 'MRS_voxel_mask_CSF';
        O_CSFvox.dim = voxmaskvol.dim;
        O_CSFvox.dt = voxmaskvol.dt;
        O_CSFvox.mat = voxmaskvol.mat;
        CSF_voxmask_vol = CSFvol.private.dat(:,:,:) .* voxmaskvol.private.dat(:,:,:);
        O_CSFvox = spm_write_vol(O_CSFvox, CSF_voxmask_vol);
        
        % 3. Calculate a CSF-corrected i.u. value and output it to the structure
        
        GMsum  = sum(sum(sum(O_GMvox.private.dat(:,:,:))));
        WMsum  = sum(sum(sum(O_WMvox.private.dat(:,:,:))));
        CSFsum = sum(sum(sum(O_CSFvox.private.dat(:,:,:))));
        
        fGM  = GMsum / (GMsum + WMsum + CSFsum);
        fWM  = WMsum / (GMsum + WMsum + CSFsum);
        fCSF = CSFsum / (GMsum + WMsum + CSFsum);
        
        MRS_struct.out.(vox{kk}).tissue.fGM(ii)  = fGM;
        MRS_struct.out.(vox{kk}).tissue.fWM(ii)  = fWM;
        MRS_struct.out.(vox{kk}).tissue.fCSF(ii) = fCSF;
        
        % Correction of institutional units only feasible if water-scaling
        % is performed, skip otherwise
        if strcmp(MRS_struct.p.reference,'H2O')
            target = [MRS_struct.p.target, {'Cr'}, {'Cho'}, {'NAA'}]; % Add Cr, Cho, and NAA
            for jj = 1:length(target)
                if strcmp(target{jj},'GABAGlx')
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
        subplot(2,3,4:6);
        axis off;
        
        text_pos = 1;
        
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{1,ii});
        end
        fname = [tmp2 tmp3];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, text_pos-0.12, 'Filename: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, [' ' fname], 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        [~,tmp2,tmp3] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp2 tmp3];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        text(0.5, text_pos-0.24, 'Anatomical image: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, [' ' T1image], 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13, 'Interpreter', 'none');
        
        % Print correction of institutional units only feasible if water-scaling is performed, skip otherwise
        text_pos = text_pos-0.24;
        if strcmp(MRS_struct.p.reference, 'H2O')
            target = MRS_struct.p.target;
            tmp = strcmp(target, 'GABAGlx');
            if any(tmp)
                if MRS_struct.p.HERMES
                    target = {'GABA','Glx',target{~tmp}};
                else
                    target = {'GABA','Glx'};
                end
            end
            for jj = 1:length(target)
                text_pos = text_pos - 0.12;
                switch target{jj}
                    case 'GABA'
                        tmp1 = 'GABA+/Water (CSF-corrected): ';
                        tmp2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU_CSFcorr(ii));
                    case 'Lac'
                        tmp1 = 'Lac+MM/Water (CSF-corrected): ';
                        tmp2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).Lac.ConcIU_CSFcorr(ii));
                    case {'Glx','GSH','EtOH'}
                        tmp1 = [target{jj} '/Water (CSF-corrected): '];
                        tmp2 = sprintf(' %.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_CSFcorr(ii));
                end
                text(0.5, text_pos, tmp1, 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
                text(0.5, text_pos, tmp2, 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
            end
        end
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fGM(ii));
        text(0.5, text_pos-0.12, 'GM voxel fraction: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.12, tmp, 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fWM(ii));
        text(0.5, text_pos-0.24, 'WM voxel fraction: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.24, tmp, 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        tmp = sprintf(' %.2f', MRS_struct.out.(vox{kk}).tissue.fCSF(ii));
        text(0.5, text_pos-0.36, 'CSF voxel fraction: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.36, tmp, 'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        text(0.5, text_pos-0.48, 'SegmentVer: ', 'FontName', 'Arial', 'HorizontalAlignment','right', 'VerticalAlignment', 'top', 'FontSize', 13);
        text(0.5, text_pos-0.48, [' ' MRS_struct.version.segment],  'FontName', 'Arial', 'VerticalAlignment', 'top', 'FontSize', 13);
        
        % Voxel segmentation
        if isfield(MRS_struct.p,'TablePosition')
            voxoff = MRS_struct.p.voxoff(ii,:) + MRS_struct.p.TablePosition(ii,:);
        else
            voxoff = MRS_struct.p.voxoff(ii,:);
        end
        
        % Plot segmented voxel as a montage
        MRS_struct.mask.(vox{kk}).img_montage{ii} = PlotSegmentedVoxels(struc, voxoff, voxmaskvol, O_GMvox, O_WMvox, O_CSFvox);

        if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii});
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
        if exist(MRS_struct.out.(vox{kk}).CSVname, 'file')
            fprintf('\nUpdating results in %s\n', [MRS_struct.out.(vox{kk}).CSVname '...']);
        else
            fprintf('\nExporting results to %s\n', [MRS_struct.out.(vox{kk}).CSVname '...']);
        end
    end
    
end

warning('on'); % turn warnings back on

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end


function img_montage = PlotSegmentedVoxels(struc, voxoff, voxmaskvol, O_GMvox, O_WMvox, O_CSFvox)

img_t      = flipud(voxel2world_space(spm_vol(struc), voxoff));
mask_t     = flipud(voxel2world_space(voxmaskvol, voxoff));
mask_t_GM  = flipud(voxel2world_space(O_GMvox, voxoff));
mask_t_WM  = flipud(voxel2world_space(O_WMvox, voxoff));
mask_t_CSF = flipud(voxel2world_space(O_CSFvox, voxoff));

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

ha = subplot(2,3,1:3);
imagesc(img_montage);
axis equal tight off;
text(floor(size(mask_t,2)/2), 20, 'Voxel', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'GM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(2*floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'WM', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
text(3*floor(size(mask_t,2)) + floor(size(mask_t,2)/2), 20, 'CSF', 'Color', [1 1 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
set(ha, 'pos', [0 0.17 1 1]);



