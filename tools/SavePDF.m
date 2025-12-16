function run_count = SavePDF(h, MRS_struct, ii, jj, kk, vox, module, run_count, target)

if ~isMATLABReleaseOlderThan("R2025a") && MRS_struct.p.append
    font_size_adj = 2.75;
else
    font_size_adj = 0;
end

% Gannet logo
axes('Position', [0.8825, 0.04, 0.125, 0.125], 'Units', 'normalized');
Gannet_logo = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'misc', 'Gannet3_logo.png');
I = imread(Gannet_logo, 'BackgroundColor', 'none');
imshow(I);
axis off image;

% Gannet version
d.left   = 0;
d.bottom = 0.02;
d.width  = 1;
d.height = 0.02;
axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
text(0.9925, 0, MRS_struct.info.version.Gannet, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
axis off;

% Gannet documentation
axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
str = 'For complete documentation, please visit: https://markmikkelsen.github.io/Gannet-docs';
text(0.5, 0, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11 - font_size_adj, 'HorizontalAlignment', 'center');
axis off square;

% Batch number and output time
d.bottom = 0.98;
axes('Position', [d.left d.bottom d.width d.height], 'Units', 'normalized');
if strcmp(module,'GannetFit') && MRS_struct.p.HERMES && MRS_struct.p.append && ~isempty(fileparts(which('export_fig')))
    if strcmp(target,'GABAGlx')
        target = 'GABA+Glx';
    end
    text(0.0075, 0, ['Batch file: ' num2str(ii) ' of ' num2str(MRS_struct.p.numScans) ' (' target ')'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11 - font_size_adj, 'HorizontalAlignment', 'left');
else
    text(0.0075, 0, ['Batch file: ' num2str(ii) ' of ' num2str(MRS_struct.p.numScans)], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11 - font_size_adj, 'HorizontalAlignment', 'left');
end
text(0.9925, 0, char(datetime('now','Format','dd-MMM-y HH:mm:ss')), 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 11 - font_size_adj, 'HorizontalAlignment', 'right');
axis off;

if any(strcmp(listfonts, 'Arial'))
    set(findall(h, '-property', 'FontName'), 'FontName', 'Arial');
end
set(findall(h, '-property', 'XColor'), 'XColor', [0 0 0]);
set(findall(h, '-property', 'YColor'), 'YColor', [0 0 0]);

% If export_fig is installed, export PDF using it
if MRS_struct.p.append && ~isempty(fileparts(which('export_fig')))

    scr_sz = get(0,'ScreenSize');
    if ispc
        px_sz = 96;
    elseif ismac || isunix
        px_sz = 72;
    end
    fig_w = 11*px_sz;
    fig_h = 8.5*px_sz;
    set(h, 'Units', 'Pixels', 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);

    % Create output folder
    if ~MRS_struct.p.bids
        out_dir = fullfile(pwd, 'Gannet_output');
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    else % BIDSify
        out_dir = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', 'pdfs');
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end

    pdf_name = fullfile(out_dir, [module '.pdf']);
    if exist(pdf_name, 'file') && (ii + jj) == 2
        run_count = 1;
        pdf_name  = fullfile(out_dir, [module num2str(run_count) '.pdf']);
        while 1
            if exist(pdf_name, 'file')
                run_count = run_count + 1;
                pdf_name  = fullfile(out_dir, [module num2str(run_count) '.pdf']);
            else
                break
            end
        end
    elseif (ii + jj) > 2 && run_count > 0
        pdf_name = fullfile(out_dir, [module num2str(run_count) '.pdf']);
    end

    export_fig(pdf_name, '-pdf', '-painters', '-append', '-nocrop', '-nofontswap', '-silent', h);

else

    if MRS_struct.p.append && isempty(fileparts(which('export_fig'))) && ii == 1
        warning(['Could not find the function ''export_fig.m''. ', ...
                 'Cannot append PDFs. ', ...
                 'Please ensure that you have added the export_fig/ ', ...
                 'folder in the main Gannet folder to your MATLAB ', ...
                 'search path. PDFs will be saved separately.']);
    end

    set(h, 'PaperUnits', 'inches', 'PaperSize', [11 8.5], 'PaperPosition', [0 0 11 8.5]);

    % Create output folder
    if ~MRS_struct.p.bids
        out_dir = fullfile(pwd, [module '_output']);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    else % BIDSify
        out_dir = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', 'pdfs', module);
        if ~exist(out_dir, 'dir')
            mkdir(out_dir);
        end
    end

    % For Philips .data
    if strcmp(MRS_struct.p.vendor, 'Philips_data')
        fullpath = MRS_struct.metabfile{1,ii};
        fullpath = regexprep(fullpath, '.data', '_data');
        fullpath = regexprep(fullpath, '\', '_');
        fullpath = regexprep(fullpath, '/', '_');
    end

    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
        [~, metabfile_nopath] = fileparts(MRS_struct.metabfile{1,ii*2-1});
    else
        [~, metabfile_nopath, ext] = fileparts(MRS_struct.metabfile{1,ii});
        if strcmpi(ext, '.gz')
            metabfile_nopath(end-3:end) = [];
        end
    end

    module2 = lower(module);
    module2(1:6) = [];

    if strcmp(module2, 'fit') && MRS_struct.p.HERMES
        module2 = [MRS_struct.p.target{jj} '_' module2];
    end

    if strcmp(MRS_struct.p.vendor, 'Philips_data')
        if isfield(MRS_struct.p, 'trimmed_avgs')
            pdf_name = fullfile(out_dir, [fullpath '_' vox{kk} '_' module2 '_' num2str(MRS_struct.p.Navg(ii)) '_avgs.pdf']);
        else
            pdf_name = fullfile(out_dir, [fullpath '_' vox{kk} '_' module2 '.pdf']);
        end
    else
        if isfield(MRS_struct.p, 'trimmed_avgs')
            pdf_name = fullfile(out_dir, [metabfile_nopath '_' vox{kk} '_' module2 '_' num2str(MRS_struct.p.Navg(ii)) '_avgs.pdf']);
        else
            pdf_name = fullfile(out_dir, [metabfile_nopath '_' vox{kk} '_' module2 '.pdf']);
        end
    end

    saveas(h, pdf_name);

end



