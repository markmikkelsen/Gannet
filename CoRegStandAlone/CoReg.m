function MRS_struct = CoReg(MRS_struct, struc)

% Coregistration of MRS voxel volumes to imaging datasets, based on headers.

loadFile = which('GannetCoRegister');
fileID = fopen(loadFile, 'rt');
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.coreg = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
MRS_struct.version.coreg = out.version;

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

if MRS_struct.ii ~= length(struc)
    error('The number of NIfTI files does not match the number of MRS files processed by CoRegStandAlone.');
end

numscans = numel(MRS_struct.metabfile);
vox = MRS_struct.p.vox(1);

for ii = 1:numscans
    
    % Loop over voxels if PRIAM
    for kk = 1:length(vox)
        
        %Ultimately this switch will not be necessary...
        switch MRS_struct.p.vendor
            
            case 'Philips'
                sparname = [MRS_struct.metabfile{ii}(1:(end-4)) MRS_struct.p.spar_string];
                MRS_struct = GannetMask_Philips(sparname, struc{ii}, MRS_struct, ii, vox, kk);
                
            case 'Philips_data'
                if exist(MRS_struct.metabfile_sdat,'file')
                    MRS_struct.p.vendor = 'Philips';
                    MRS_struct.metabfile_data = MRS_struct.metabfile;
                    MRS_struct.metabfile = MRS_struct.metabfile_sdat;
                    MRS_struct = GannetCoRegister(MRS_struct, struc);
                    MRS_struct.metabfile = MRS_struct.metabfile_data;
                    MRS_struct.p.vendor = 'Philips_data';
                else
                    error([MRS_struct.p.vendor ' format does not include voxel location information in the header. See notes in GannetCoRegister.']);
                    %If this comes up, once GannetLoad has been read:
                    %1. Switch vendor to Philips
                    %       MRS_struct.p.vendor = 'Philips';
                    %2. Copy .data filenames.
                    %       MRS_struct.metabfile_data = MRS_struct.metabfile;
                    %3. Replace the list with the corrsponding SDAT files (in correct order)
                    %        MRS_struct.metabfile = {'SDATfile1.sdat' 'SDATfile2.SDAT'};
                    %4. Rerun GannetCoRegister
                    %
                    %5.  Copy .sdat filenames and replace .data ones. Tidy up.
                    %       MRS_struct.metabfile_sdat = MRS_struct.metabfile;
                    %       MRS_struct.metabfile = MRS_struct.metabfile_data;
                    %       MRS_struct.p.vendor = 'Philips_data'
                end
                
            case 'Siemens_rda'
                MRS_struct = GannetMask_SiemensRDA(MRS_struct.metabfile{ii}, struc{ii}, MRS_struct, ii, vox, kk);
                
            case {'Siemens_twix', 'Siemens_dicom', 'dicom'}
                MRS_struct = GannetMask_SiemensTWIX(MRS_struct.metabfile{ii}, struc{ii}, MRS_struct, ii, vox, kk);
                
            case 'GE'
                [~,~,ext] = fileparts(struc{ii});
                if strcmp(ext,'.nii')
                    MRS_struct = GannetMask_GE_nii(MRS_struct.metabfile{ii}, struc{ii}, MRS_struct, ii, vox, kk);
                else
                    MRS_struct = GannetMask_GE(MRS_struct.metabfile{ii}, struc{ii}, MRS_struct, ii, vox, kk);
                end

            case 'nifti'
                error('NIfTI not yet supported.');
                
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
        scr_sz = get(0,'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetCoRegister Output';
        set(gcf,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
        
        subplot(2,3,4:6);
        axis off;
        
        [~,tmp,tmp2] = fileparts(MRS_struct.mask.(vox{kk}).outfile{ii});
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.5, 0.75, 'Mask output: ', 'HorizontalAlignment' , 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.75, [' ' fname], 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'none');
        
        text(0.5, 0.63, 'Spatial parameters: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.63, ' [LR, AP, FH]', 'FontName', 'Arial', 'FontSize', 13);
        
        tmp = [' ' num2str(MRS_struct.p.voxdim(ii,1)) ' \times ' num2str(MRS_struct.p.voxdim(ii,2)) ' \times ' num2str(MRS_struct.p.voxdim(ii,3)) ' mm^{3}'];
        text(0.5, 0.51, 'Dimensions: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.51, tmp, 'FontName', 'Arial', 'FontSize', 13, 'Interpreter', 'tex');
        
        tmp = [' ' num2str(prod(MRS_struct.p.voxdim(ii,:))/1e3) ' mL'];
        text(0.5, 0.39, 'Volume: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.39, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        tmp = [' [' num2str(MRS_struct.p.voxoff(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxoff(ii,3), '%3.1f') '] mm'];
        text(0.5, 0.27, 'Position: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.27, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        if any(strcmp(MRS_struct.p.vendor,{'Philips','Philips_data'}))
            tmp = [' [' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        else
            tmp = [' [' num2str(MRS_struct.p.voxang(ii,1), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,2), '%3.1f') ', ' num2str(MRS_struct.p.voxang(ii,3), '%3.1f') '] deg'];
        end
        text(0.5, 0.15, 'Angulation: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.15, tmp, 'FontName', 'Arial', 'FontSize', 13);
        
        text(0.5, 0.03, 'CoRegVer: ', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 13);
        text(0.5, 0.03, [' ' MRS_struct.version.coreg], 'FontName', 'Arial', 'FontSize', 13);
        
        ha = subplot(2,3,1:3);
        
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
        
        imagesc(squeeze(MRS_struct.mask.(vox{kk}).img{ii}));
        colormap('gray');
        img = MRS_struct.mask.(vox{kk}).img{ii}(:);
        caxis([0 mean(img(img > 0.01)) + 3*std(img(img > 0.01))]);
        axis equal tight off;
        text(10, size(MRS_struct.mask.(vox{kk}).img{ii},1)/2, 'L', 'Color', [1 1 1], 'FontSize', 20);
        text(size(MRS_struct.mask.(vox{kk}).img{ii},2) - 20, size(MRS_struct.mask.(vox{kk}).img{ii},1)/2, 'R', 'Color', [1 1 1], 'FontSize', 20);
        set(ha,'pos',[0 0.15 1 1]);
        title(t, 'FontName', 'Arial', 'FontSize', 15, 'Interpreter', 'none');
        
        % Gannet logo
        Gannet_logo = fullfile(fileparts(which('GannetLoad')), 'Gannet3_logo.png');
        I = imread(Gannet_logo,'png','BackgroundColor',[1 1 1]);
        axes('Position', [0.825, 0.05, 0.125, 0.125]);
        imshow(I);
        text(0.9, 0, MRS_struct.version.Gannet, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
        axis off square;
        
        % Gannet documentation
        axes('Position', [(1-0.9)/2, 0.025, 0.9, 0.15]);
        str = 'For complete documentation, please visit: https://markmikkelsen.github.io/Gannet-docs';
        text(0.5, 0, str, 'FontName', 'Arial', 'FontSize', 11, 'HorizontalAlignment', 'center');
        axis off square;
        
        % For Philips .data
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            fullpath = MRS_struct.metabfile{ii};
            fullpath = regexprep(fullpath, '.data', '_data');
            fullpath = regexprep(fullpath, '\', '_');
            fullpath = regexprep(fullpath, '/', '_');
        end
        
        [~,metabfile_nopath] = fileparts(MRS_struct.metabfile{ii});
        
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
        set(gcf,'PaperUnits','inches');
        set(gcf,'PaperSize',[11 8.5]);
        set(gcf,'PaperPosition',[0 0 11 8.5]);
        if strcmpi(MRS_struct.p.vendor,'Philips_data')
            pdfname = fullfile(pwd, 'CoRegStandAlone_output', [fullpath '_' vox{kk} '_coreg.pdf']);
        else
            pdfname = fullfile(pwd, 'CoRegStandAlone_output', [metabfile_nopath '_' vox{kk} '_coreg.pdf']);
        end
        saveas(gcf, pdfname);
        
    end
    
end

warning('on'); % turn warnings back on

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide
    close(figTitle);
end



