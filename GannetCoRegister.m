function MRS_struct = GannetCoRegister(MRS_struct, struc)
% Co-registration of MRS voxel volumes to imaging datasets, based on headers.

if nargin < 1 || ...
    (nargin < 2 && ~(isfield(MRS_struct, 'p') && isfield(MRS_struct.p, 'bids') && MRS_struct.p.bids))
    fprintf('\n');
    error('MATLAB:minrhs', ['Not enough input arguments. ' ...
          'GannetCoRegister requires two arguments (MRS_struct, struc), ' ...
          'unless processing a BIDS dataset, in which case only MRS_struct ' ...
          'is required.']);
end

if ~isstruct(MRS_struct)
    fprintf('\n');
    error('The first input argument must be a structure, but received %s.', class(MRS_struct));
end

if nargin == 2
    if ~iscell(struc)
        fprintf('\n');
        error('The second input argument ''%s'' must be a structure.', struc);
    end
end

MRS_struct.info.datetime.coreg = datetime('now');
MRS_struct.info.version.coreg = '250912';

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

% Find MR images if processing a BIDS dataset
if MRS_struct.p.bids
    struc = cell(MRS_struct.p.numScans,1);
    for ii = 1:MRS_struct.p.numScans
        bids_file = bids.File(MRS_struct.metabfile{ii});
        if ~exist(fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path), 'dir')
            bids.util.mkdir(fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path));
        end
        metadata = bids.internal.get_metadata(bids.internal.get_meta_list(MRS_struct.metabfile{ii}));
        try
            struc{ii} = bids.internal.resolve_bids_uri(metadata.AnatomicalImage, MRS_struct.out.BIDS);
        catch
            fprintf('\n');
            error(['No valid structural images found for ''%s''.' ...
                   '\nCheck that its JSON sidecar file has an entry for ''AnatomicalImage''.'], bids_file.filename);
        end
    end
else
    struc = GetFullPath(struc);
end

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
        if ~isMATLABReleaseOlderThan("R2025a")
            h.Theme = 'light';
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

        [~,name,ext] = fileparts(MRS_struct.mask.(vox{kk}).fname{ii});
        fname = [name ext];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, 0.75, 'Mask output: ', 'Units', 'normalized', 'HorizontalAlignment' , 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.75, [' ' fname], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');

        text(0.5, 0.63, 'Spatial parameters: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.63, ' [LR, PA, SI]', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        str = [' ' num2str(MRS_struct.p.voxdim(ii,1)) ' \times ' num2str(MRS_struct.p.voxdim(ii,2)) ' \times ' num2str(MRS_struct.p.voxdim(ii,3)) ' mm^{3}'];
        text(0.5, 0.51, 'Dimensions: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.51, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'tex');

        str = [' ' num2str(prod(MRS_struct.p.voxdim(ii,:))/1e3) ' mL'];
        text(0.5, 0.39, 'Volume: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.39, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        str = [' [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'];
        text(0.5, 0.27, 'Position: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.27, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        if any(strcmp(MRS_struct.p.vendor, {'Philips', 'Philips_data'}))
            str = [' [' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        else
            str = [' [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        end
        text(0.5, 0.15, 'Angulation: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.15, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        text(0.5, 0.03, 'CoRegVer: ', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.03, [' ' MRS_struct.info.version.coreg], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 13);

        hb = subplot(2,3,1:3);

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

        imagesc(MRS_struct.mask.(vox{kk}).img{ii});
        axis equal tight off;
        text(0.01, 0.5, 'L', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
        text(0.16, 0.95, 'A', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
        text(0.32, 0.5, 'A', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
        text(0.5, 0.95, 'S', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
        text(0.825, 0.95, 'S', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
        text(0.975, 0.5, 'R', 'Color', [1 1 1], 'FontSize', 20, 'Units', 'normalized');
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



