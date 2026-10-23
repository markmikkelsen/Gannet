function PlotAvgSubpec(MRS_struct, ii, kk, vox, n_subspec, experiment)
% Plot and print the averaged ON/OFF (MEGA) or A–D subspectra (HERMES)
% This is mostly useful for debugging or assessing data quality issues

spec_plot_dir = fullfile(pwd, 'Gannet_spec_plots');
if ~exist(spec_plot_dir, 'dir')
    mkdir(spec_plot_dir);
end

if ishandle(401)
    clf(401);
end
if MRS_struct.p.hide
    h = figure('Visible', 'off');
else
    h = figure(401);
end
if ~isMATLABReleaseOlderThan("R2025a")
    h.Theme = 'light';
end

% First, plot subspectra separately
figDims.w = 0.7;
if MRS_struct.p.HERMES
    figDims.h = 0.7;
else
    figDims.h = 0.5;
end
figDims.l = (1 - figDims.w)/2;
figDims.b = (1 - figDims.h)/2;
set(h, 'Color', 'w', 'Units', 'Normalized', ...
    'OuterPosition', [figDims.l figDims.b figDims.w figDims.h]);

p = tiledlayout('flow','TileSpacing','compact','Padding','compact');
for ll = 1:n_subspec
    ax = nexttile(p);
    plot(ax, MRS_struct.spec.freq, ...
        real(MRS_struct.spec.(vox{kk}).subspec.(experiment{ll})(ii,:)));
    set(ax, 'XDir', 'reverse', 'XLim', [0.5 4.25], 'box', 'off', ...
        'TickDir', 'out');
    title((experiment{ll}), 'FontSize', 20);
    xlabel('ppm', 'FontSize', 16);
end

[~,fname] = fileparts(MRS_struct.metabfile{ii});
if MRS_struct.p.HERMES
    png_name = fullfile(spec_plot_dir, [fname '_A-D_subspec_separate.png']);
else
    png_name = fullfile(spec_plot_dir, [fname '_ON_OFF_subspec_separate.png']);
end
exportgraphics(h, png_name, 'Resolution', 300);

% Now plot them overlapped
clf(h);
figDims.w = 0.5;
figDims.h = 0.5;
figDims.l = (1 - figDims.w)/2;
figDims.b = (1 - figDims.h)/2;
set(h, 'OuterPosition', [figDims.l figDims.b figDims.w figDims.h]);

hold on;
for ll = 1:n_subspec
    plot(MRS_struct.spec.freq, ...
        real(MRS_struct.spec.(vox{kk}).subspec.(experiment{ll})(ii,:)));
    set(gca, 'XDir', 'reverse', 'XLim', [0.5 4.25], 'box', 'off', ...
        'TickDir', 'out');
    xlabel('ppm', 'FontSize', 16);
end
hold off;
legend(experiment, 'FontSize', 20, 'Box', 'off', 'Location', 'northeast');

if MRS_struct.p.HERMES
    png_name = fullfile(spec_plot_dir, [fname '_A-D_subspec_overlapped.png']);
else
    png_name = fullfile(spec_plot_dir, [fname '_ON_OFF_subspec_overlapped.png']);
end
exportgraphics(h, png_name, 'Resolution', 300);

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide
    close(h);
end

end
