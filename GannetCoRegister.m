function MRS_struct = GannetCoRegister(MRS_struct, struc)

% Co-registration of MRS voxel volumes to imaging datasets, based on headers.

if nargin < 2
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

MRS_struct.version.coreg = '230823';

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

struc = GetFullPath(struc);

if MRS_struct.p.numScans ~= length(struc)
    fprintf('\n');
    error('The number of structural image files does not match the number of MRS files processed by GannetLoad.');
end

missing = 0;
for filecheck = 1:numel(struc)
    if ~exist(struc{filecheck}, 'file')
        fprintf('\nThe file ''%s'' (#%d) is missing. Typo?\n', struc{filecheck}, filecheck);
        missing = 1;
    end
end

if missing
    fprintf('\n');
    error('Not all structural image files could be found. Please check filenames. Exiting...');
end

run_count = 0;

for ii = 1:MRS_struct.p.numScans

    [~,b,c] = fileparts(MRS_struct.metabfile{1,ii});
    [~,e,f] = fileparts(struc{ii});
    if ii == 1
        fprintf('\nCo-registering voxel from %s to %s...\n', [b c], [e f]);
    else
        fprintf('Co-registering voxel from %s to %s...\n', [b c], [e f]);
    end

    fname = MRS_struct.metabfile{1,ii};
    if strcmpi(f, '.gz')
        fprintf('Uncompressing %s...\n', struc{ii});
        struc(ii) = gunzip(struc{ii});
    end

    % Loop over voxels if PRIAM
    for kk = 1:length(vox)

        switch MRS_struct.p.vendor

            case 'GE'
                [~,~,ext] = fileparts(struc{ii});
                if strcmpi(ext, '.nii')
                    MRS_struct = GannetMask_GE_nii(fname, struc{ii}, MRS_struct, ii, vox, kk);
                else
                    MRS_struct = GannetMask_GE(fname, struc{ii}, MRS_struct, ii, vox, kk);
                end

            case 'NIfTI'
                MRS_struct = GannetMask_NIfTI(fname, struc{ii}, MRS_struct, ii, vox, kk);

            case 'Philips'
                MRS_struct = GannetMask_Philips(fname, struc{ii}, MRS_struct, ii, vox, kk);

            case 'Philips_data'
                if exist(MRS_struct.metabfile_sdat, 'file')
                    MRS_struct.p.vendor = 'Philips';
                    MRS_struct.metabfile_data = MRS_struct.metabfile;
                    MRS_struct.metabfile = MRS_struct.metabfile_sdat;
                    MRS_struct = GannetCoRegister(MRS_struct, struc);
                    MRS_struct.metabfile = MRS_struct.metabfile_data;
                    MRS_struct.p.vendor = 'Philips_data';
                else
                    error('%s format does not include voxel location information in the header. See notes in GannetCoRegister.', MRS_struct.p.vendor);
                    % If this comes up, once GannetLoad has been read:
                    % 1. Switch vendor to Philips
                    %        MRS_struct.p.vendor = 'Philips';
                    % 2. Copy .data filenames.
                    %        MRS_struct.metabfile_data = MRS_struct.metabfile;
                    % 3. Replace the list with the corrsponding SDAT files (in correct order)
                    %         MRS_struct.metabfile = {'SDATfile1.sdat' 'SDATfile2.SDAT'};
                    % 4. Rerun GannetCoRegister
                    % 5. Copy .sdat filenames and replace .data ones. Tidy up.
                    %        MRS_struct.metabfile_sdat = MRS_struct.metabfile;
                    %        MRS_struct.metabfile = MRS_struct.metabfile_data;
                    %        MRS_struct.p.vendor = 'Philips_data'
                end

            case 'Siemens_rda'
                fname      = MRS_struct.metabfile{1,ii*2-1};
                MRS_struct = GannetMask_SiemensRDA(fname, struc{ii}, MRS_struct, ii, vox, kk);

            case {'Siemens_dicom', 'Siemens_twix', 'DICOM'}
                MRS_struct = GannetMask_SiemensTWIX(fname, struc{ii}, MRS_struct, ii, vox, kk);

        end

        % Build output figure
        if ishandle(103)
            clf(103);
        end
        if MRS_struct.p.hide
            h = figure('Visible', 'off');
        else
            h = figure(103);
        end
        % Open figure in center of screen
        scr_sz = get(0,'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetCoRegister Output';
        set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');

        ha  = subplot(2,3,4:6);
        pos = get(ha, 'Position');
        set(ha, 'Position', [0 pos(2) 1 pos(4)]);
        axis off;

        [~,tmp,tmp2] = fileparts(MRS_struct.mask.(vox{kk}).outfile{ii});
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, 0.75, 'Mask output: ', 'Units', 'normalized', 'HorizontalAlignment' , 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.75, [' ' fname], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');

        text(0.5, 0.63, 'Spatial parameters: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.63, ' [LR, AP, FH]', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        tmp = [' ' num2str(MRS_struct.p.voxdim(ii,1)) ' \times ' num2str(MRS_struct.p.voxdim(ii,2)) ' \times ' num2str(MRS_struct.p.voxdim(ii,3)) ' mm^{3}'];
        text(0.5, 0.51, 'Dimensions: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.51, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'tex');

        tmp = [' ' num2str(prod(MRS_struct.p.voxdim(ii,:))/1e3) ' mL'];
        text(0.5, 0.39, 'Volume: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.39, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        tmp = [' [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'];
        text(0.5, 0.27, 'Position: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.27, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        if any(strcmp(MRS_struct.p.vendor, {'Philips', 'Philips_data'}))
            tmp = [' [' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        else
            tmp = [' [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        end
        text(0.5, 0.15, 'Angulation: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.15, tmp, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        text(0.5, 0.03, 'CoRegVer: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.03, [' ' MRS_struct.version.coreg], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        hb = subplot(2,3,1:3);

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

        imagesc(MRS_struct.mask.(vox{kk}).img{ii});
        axis equal tight off;
        text(10, size(MRS_struct.mask.(vox{kk}).img{ii},1)/2, 'L', 'Color', [1 1 1], 'FontSize', 20);
        text(size(MRS_struct.mask.(vox{kk}).img{ii},2) - 20, size(MRS_struct.mask.(vox{kk}).img{ii},1)/2, 'R', 'Color', [1 1 1], 'FontSize', 20);
        set(hb,'Position',[0 0.15 1 1]);
        title(t, 'FontName', 'Arial', 'FontSize', 15, 'Interpreter', 'none');

        % Save output as PDF
        run_count = SavePDF(h, MRS_struct, ii, 1, kk, vox, mfilename, run_count);

    end

end

% Save MRS_struct as mat file
if MRS_struct.p.mat
    mat_name = fullfile(pwd, ['MRS_struct_' vox{kk} '.mat']);
    if exist(mat_name, 'file')
        fprintf('\nUpdating results in %s\n', ['MRS_struct_' vox{kk} '.mat...']);
    else
        fprintf('\nSaving results to %s\n', ['MRS_struct_' vox{kk} '.mat...']);
    end
    save(mat_name, 'MRS_struct', '-v7.3');
end

warning('on'); % turn warnings back on

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end



