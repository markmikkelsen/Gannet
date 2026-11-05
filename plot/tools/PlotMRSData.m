function varargout = PlotMRSData(MRS_struct, data, opts)
% Quick plot of FIDs or spectra in MRS_struct in the frequency domain. This
% is a debugging/inspection shortcut that can be called at any point in the
% Gannet pipeline
%
% Time-domain data (FIDs) are line-broadened (optional), zero-filled and
% Fourier-transformed before plotting. The ppm axis is taken from
% MRS_struct.spec.freq if it exists and matches the data; otherwise it is
% created on the fly using the same convention as GannetLoad
%
% Figure visibility is controlled by MRS_struct.p.hide (set in
% GannetPreInitialise): if hide = 1, the figure is made invisible (and
% closed after printing, unless the figure handle is requested). If
% MRS_struct.p.hide does not exist (e.g., when used standalone), the figure
% is visible
%
% Usage:
%   PlotMRSData(MRS_struct, 'fids.data', 'ii', ii)
%   PlotMRSData(MRS_struct, AllFramesFT, 'ii', ii, 'Average', true)
%   PlotMRSData(MRS_struct, 'spec.{vox}.subspec.{subspec}', 'ii', ii, 'kk', kk, 'll', 2)
%   PlotMRSData(MRS_struct, {'spec.{vox}.subspec.ON', 'spec.{vox}.subspec.OFF'}, ...
%               'ii', ii, 'kk', kk, 'Print', false)
%   PlotMRSData(MRS_struct, 'spec.{vox}.{target}.diff', 'ii', ii, 'jj', jj, 'kk', kk)
%   h = PlotMRSData(...)
%
% Inputs:
%   MRS_struct - Gannet data structure
%   data       - Data to plot. One of:
%                  - a numeric vector or matrix (e.g., a local variable such
%                    as AllFramesFT or AllFramesFTrealign)
%                  - a char/string field path relative to MRS_struct, e.g.,
%                    'fids.data' or 'spec.vox1.GABA.diff' (a leading
%                    'MRS_struct.' is allowed). The placeholders {vox},
%                    {target} and {subspec} are replaced with
%                    MRS_struct.p.vox{kk}, MRS_struct.p.target{jj} and the
%                    ll-th subspectrum name ({'ON','OFF'} for MEGA,
%                    {'A','B','C','D'} for HERMES/HERCULES)
%                  - a cell array of any of the above; these are overlaid
%                    (or plotted in separate tiles; see Layout)
%   Optional name-value arguments:
%     ii         - dataset index (default: 1); sets the acquisition
%                  parameters (sw, LarmorFreq, ZeroFillTo, ...) and the
%                  row selected from spectra stored as (ii,:) arrays
%     jj         - target index for the {target} placeholder (default: 1)
%     kk         - voxel index for the {vox} placeholder (default: 1)
%     ll         - subspectrum index for the {subspec} placeholder
%                  (default: 1)
%     Domain     - domain of the input data: 'auto' (default), 'time' or
%                  'freq'. In 'auto', paths beginning with 'fids' are time
%                  domain and paths beginning with 'spec' are frequency
%                  domain; for numeric input, data whose magnitude peaks
%                  in the first few points are treated as FIDs
%     SpecDim    - dimension of the data that runs along time/frequency
%                  (1 or 2; default: determined from the data size)
%     Transients - indices of the transients to plot, i.e., the columns of
%                  (points x transients) data or the rows of (ii,:) arrays
%                  (default: all transients for (points x transients) data;
%                  row ii for (ii,:) arrays)
%     Average    - average the selected transients before plotting
%                  (default: false)
%     ZeroFill   - zero-fill FIDs before the FFT (default: true); uses
%                  MRS_struct.p.ZeroFillTo(ii) if it exists, otherwise the
%                  same rule as GannetLoad
%     LB         - exponential line broadening applied to FIDs, in Hz
%                  (default: 0)
%     Phase0     - zero-order phase applied before plotting, in degrees
%                  (default: 0)
%     Part       - 'real' (default), 'imag' or 'abs'
%     Freq       - user-supplied ppm axis (overrides everything else)
%     SW         - spectral width (Hz) used to build the ppm axis (default:
%                  MRS_struct.p.sw(ii), or MRS_struct.p.sw_water(ii) for
%                  water data if it exists)
%     F0         - ppm value at the center of the spectrum (default: 4.68,
%                  or 4.8 for phantom data)
%     XLim       - ppm limits (default: [0.5 4.7], or [3.7 5.7] for water
%                  data; [] shows the full ppm range)
%     YLim       - y-axis limits (default: for non-water data, the upper
%                  limit is 15% above the maximum signal between 0.5 and
%                  4.1 ppm, so that residual water does not dominate the
%                  plot, and the lower limit is the minimum signal in that
%                  range minus 5% of the signal range; for water data, the
%                  data within XLim plus 5% padding)
%     Layout     - 'overlay' (default) or 'tiled' (one tile per input)
%     Legend     - legend labels for multiple inputs, also used as the tile
%                  titles in the 'tiled' layout (default: field paths)
%     Title      - figure title (default: file name and field path; ''
%                  for no title)
%     LineWidth  - line width (default: 1)
%     Axes       - existing axes to plot into (hold on); no new figure is
%                  made
%     FigNum     - figure number to (re)use (default: a new figure);
%                  ignored for hidden figures, which are always new
%     Print      - print the figure as a PNG (default: true)
%     OutDir     - folder to save the PNG in; created if it does not exist
%                  (default: a subject-specific folder on the Desktop,
%                  named after the data file, e.g., ~/Desktop/<file>/)
%     FileName   - file name of the PNG (default: built from the data file
%                  name, field path and a time stamp)
%     Resolution - PNG resolution in dpi (default: 300)
%     Verbose    - print the path of the saved PNG to the command window
%                  (default: false)
%
% Output:
%   h - figure handle (optional)
%
% Author: Mark Mikkelsen, Ph.D.
% Created using Claude Opus 5.5 (2026-11-04); human-inspected by the author

arguments
    MRS_struct struct
    data
    opts.ii (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.jj (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.kk (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.ll (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Domain (1,:) char {mustBeMember(opts.Domain, {'auto', 'time', 'freq'})} = 'auto'
    opts.SpecDim double {mustBeScalarOrEmpty} = []
    opts.Transients double {mustBeInteger, mustBePositive} = []
    opts.Average (1,1) logical = false
    opts.ZeroFill (1,1) logical = true
    opts.LB (1,1) double {mustBeNonnegative} = 0
    opts.Phase0 (1,1) double = 0
    opts.Part (1,:) char {mustBeMember(opts.Part, {'real', 'imag', 'abs'})} = 'real'
    opts.Freq double = []
    opts.SW double {mustBeScalarOrEmpty} = []
    opts.F0 double {mustBeScalarOrEmpty} = []
    opts.XLim double = NaN
    opts.YLim double = []
    opts.Layout (1,:) char {mustBeMember(opts.Layout, {'overlay', 'tiled'})} = 'overlay'
    opts.Legend = {}
    opts.Title = []
    opts.LineWidth (1,1) double {mustBePositive} = 1
    opts.Axes = []
    opts.FigNum double {mustBeScalarOrEmpty} = []
    opts.Print (1,1) logical = true
    opts.OutDir {mustBeTextScalar} = ''
    opts.FileName {mustBeTextScalar} = ''
    opts.Resolution (1,1) double {mustBePositive} = 300
    opts.Verbose (1,1) logical = false
end

ii = opts.ii;

% Name of the data file (for the title and the PNG file name)
[~, fname] = GetSpecPlotDir(MRS_struct, ii, false);


%%% 1. Resolve, transform and select the data %%%

if ~iscell(data)
    data = {data};
end
nInputs = numel(data);
specs   = cell(1, nInputs);
freqs   = cell(1, nInputs);
labels  = cell(1, nInputs);
isWater = false(1, nInputs);

for mm = 1:nInputs

    [d, labels{mm}, isPath] = ResolveData(MRS_struct, data{mm}, opts, mm);
    isWater(mm) = isPath && contains(labels{mm}, 'water', 'IgnoreCase', true);

    if ~isnumeric(d) || isempty(d)
        error('PlotMRSData:badData', 'Input %d (%s) is not a non-empty numeric array.', mm, labels{mm});
    end
    d = double(squeeze(d));
    if ~ismatrix(d)
        error('PlotMRSData:badData', ['Input %d (%s) has more than two non-singleton dimensions. ' ...
            'Index it down to a vector or matrix first.'], mm, labels{mm});
    end

    % Orient the data so that time/frequency runs down the columns
    [d, rowsAreDatasets] = OrientData(MRS_struct, d, ii, opts.SpecDim);

    % Determine the domain
    domain = opts.Domain;
    if strcmp(domain, 'auto')
        domain = GuessDomain(d, labels{mm}, isPath);
    end

    sw = GetSW(MRS_struct, ii, opts.SW, isWater(mm));

    % Line-broadening, zero-filling, and FFT of time-domain data
    if strcmp(domain, 'time')
        N = size(d,1);
        if opts.LB > 0
            time = (1:N)' / sw;
            d = d .* exp(-time * opts.LB * pi);
        end
        if opts.ZeroFill
            nfft = max(GetZeroFillTo(MRS_struct, ii, sw), N);
        else
            nfft = N;
        end
        d = fftshift(fft(d, nfft, 1), 1);
    end

    % Select transients
    if isempty(opts.Transients)
        if rowsAreDatasets && size(d,2) > 1 && size(d,2) >= ii
            d = d(:,ii);
        end
    else
        if any(opts.Transients > size(d,2))
            error('PlotMRSData:badTransients', 'Transients exceed the number of transients (%d) in input %d (%s).', ...
                size(d,2), mm, labels{mm});
        end
        d = d(:,opts.Transients);
    end
    if opts.Average
        d = mean(d, 2);
    end

    % Phase and select the part to plot
    if opts.Phase0 ~= 0
        d = d * exp(1i * opts.Phase0 * pi/180);
    end
    switch opts.Part
        case 'real'
            d = real(d);
        case 'imag'
            d = imag(d);
        case 'abs'
            d = abs(d);
    end

    specs{mm} = d;
    freqs{mm} = GetFreqAxis(MRS_struct, ii, size(d,1), sw, opts);

end


%%% 2. Plot %%%

% Axis limits
if isscalar(opts.XLim) && isnan(opts.XLim)
    if all(isWater)
        xLim = [3.7 5.7];
    else
        xLim = [0.5 4.7];
    end
else
    xLim = opts.XLim;
end

% Figure and axes
createdFig = isempty(opts.Axes);
if createdFig
    % Visible unless MRS_struct.p.hide exists and is set
    vis = ~(isfield(MRS_struct, 'p') && isfield(MRS_struct.p, 'hide') && MRS_struct.p.hide);
    if ~vis
        % New hidden figure (figure(FigNum) would make it visible)
        h = figure('Visible', 'off');
    elseif isempty(opts.FigNum)
        h = figure;
    else
        h = figure(opts.FigNum);
        clf(h);
    end
    if ~isMATLABReleaseOlderThan("R2025a")
        h.Theme = 'light';
    end
    set(h, 'Color', 'w');
else
    vis = true;
    ax  = opts.Axes;
    h   = ancestor(ax, 'figure');
end

tiled = strcmp(opts.Layout, 'tiled') && createdFig && nInputs > 1;
if tiled
    tl = tiledlayout(h, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
elseif createdFig
    ax = axes(h);
end

if isempty(opts.Legend)
    legLabels = labels;
else
    legLabels = cellstr(opts.Legend);
end

colors   = get(groot, 'defaultAxesColorOrder');
lineHdls = gobjects(1, nInputs);
for mm = 1:nInputs
    if tiled
        ax = nexttile(tl);
    end
    hold(ax, 'on');
    if nInputs > 1 && ~tiled
        % One color per input when overlaying several inputs
        c = colors(mod(mm-1, size(colors,1)) + 1, :);
        hl = plot(ax, freqs{mm}, specs{mm}, 'Color', c, 'LineWidth', opts.LineWidth);
    else
        hl = plot(ax, freqs{mm}, specs{mm}, 'LineWidth', opts.LineWidth);
    end
    lineHdls(mm) = hl(1);
    if tiled
        FormatAxes(ax, freqs(mm), specs(mm), xLim, opts.YLim, isWater(mm));
        title(ax, legLabels{mm}, 'FontSize', 14, 'Interpreter', 'none');
        hold(ax, 'off');
    end
end

if ~tiled
    FormatAxes(ax, freqs, specs, xLim, opts.YLim, all(isWater));
    if createdFig
        hold(ax, 'off');
    end
    if nInputs > 1 || ~isempty(opts.Legend)
        legend(lineHdls, legLabels, 'FontSize', 14, 'Box', 'off', ...
            'Location', 'northeast', 'Interpreter', 'none');
    end
end

% Title
if isnumeric(opts.Title) % default
    if nInputs == 1
        titleStr = {fname, labels{1}};
        titleStr = titleStr(~cellfun(@isempty, titleStr));
    else
        titleStr = fname;
    end
else
    titleStr = opts.Title;
end
if ~isempty(titleStr)
    if tiled
        title(tl, titleStr, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');
    else
        title(ax, titleStr, 'FontSize', 16, 'Interpreter', 'none');
    end
end
drawnow;


%%% 3. Print as PNG (to a subject-specific folder on the Desktop by default) %%%

if opts.Print
    if isempty(opts.OutDir)
        if ispc
            homeDir = getenv('USERPROFILE');
        else
            homeDir = getenv('HOME');
        end
        desktopDir = fullfile(homeDir, 'Desktop');
        if ~isfolder(desktopDir)
            warning('PlotMRSData:noDesktop', 'Desktop folder not found; saving to the current folder instead.');
            desktopDir = pwd;
        end
        outDir = fullfile(desktopDir, fname); % fname is '' if the data file is unknown
    else
        outDir = char(opts.OutDir);
    end
    if ~isfolder(outDir)
        mkdir(outDir);
    end
    if isempty(opts.FileName)
        if isempty(fname)
            base = 'PlotMRSData';
        else
            base = fname;
        end
        if nInputs == 1
            tag = ['_' regexprep(labels{1}, '[^\w]', '_')];
        else
            tag = '';
        end
        pngName = sprintf('%s%s_%s.png', base, tag, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    else
        [~, n] = fileparts(char(opts.FileName));
        pngName = [n '.png'];
    end
    pngName = fullfile(outDir, pngName);
    exportgraphics(h, pngName, 'Resolution', opts.Resolution);
    if opts.Verbose
        fprintf('Plot saved to: %s\n', pngName);
    end
end

% Close hidden figures (as elsewhere in Gannet) unless the handle is requested
if createdFig && ~vis && nargout == 0
    close(h);
end

if nargout > 0
    varargout{1} = h;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d, label, isPath] = ResolveData(MRS_struct, in, opts, mm)
% Return the numeric data and a label for one input (numeric array or field
% path relative to MRS_struct)

if isnumeric(in)
    d      = in;
    label  = sprintf('data %d', mm);
    isPath = false;
    return
end

if ~(ischar(in) || (isstring(in) && isscalar(in)))
    error('PlotMRSData:badInput', 'Input %d must be a numeric array or a field path (char/string).', mm);
end
isPath = true;

fieldPath = regexprep(char(in), '^MRS_struct\.', '');

% Replace placeholders
if contains(fieldPath, '{vox}')
    fieldPath = strrep(fieldPath, '{vox}', GetCellEntry(MRS_struct, 'vox', opts.kk, 'kk'));
end
if contains(fieldPath, '{target}')
    fieldPath = strrep(fieldPath, '{target}', GetCellEntry(MRS_struct, 'target', opts.jj, 'jj'));
end
if contains(fieldPath, '{subspec}')
    names = SubspecNames(MRS_struct);
    if opts.ll > numel(names)
        error('PlotMRSData:badIndex', 'll = %d exceeds the number of subspectra (%d).', opts.ll, numel(names));
    end
    fieldPath = strrep(fieldPath, '{subspec}', names{opts.ll});
end

% Walk the field path
parts = strsplit(fieldPath, '.');
d = MRS_struct;
for q = 1:numel(parts)
    if ~isstruct(d) || ~isfield(d, parts{q})
        error('PlotMRSData:badPath', 'Field ''%s'' not found in MRS_struct.%s', ...
            parts{q}, strjoin(parts(1:q-1), '.'));
    end
    d = d.(parts{q});
end
label = fieldPath;

end


function out = GetCellEntry(MRS_struct, fieldName, idx, idxName)
% Return MRS_struct.p.(fieldName){idx} as a char

if ~isfield(MRS_struct, 'p') || ~isfield(MRS_struct.p, fieldName)
    error('PlotMRSData:missingField', 'MRS_struct.p.%s does not exist.', fieldName);
end
list = cellstr(MRS_struct.p.(fieldName));
if idx > numel(list)
    error('PlotMRSData:badIndex', '%s = %d exceeds the number of entries (%d) in MRS_struct.p.%s.', ...
        idxName, idx, numel(list), fieldName);
end
out = list{idx};

end


function names = SubspecNames(MRS_struct)
% Subspectrum names, ordered as in SignalAveraging

if isfield(MRS_struct.p, 'HERMES') && MRS_struct.p.HERMES
    names = {'A', 'B', 'C', 'D'};
else
    names = {'ON', 'OFF'};
    if isfield(MRS_struct, 'fids') && isfield(MRS_struct.fids, 'ON_OFF') ...
            && ~MRS_struct.fids.ON_OFF(1) % OFF is first
        names = flip(names);
    end
end

end


function [d, rowsAreDatasets] = OrientData(MRS_struct, d, ii, specDim)
% Orient the data as (points x transients). Gannet stores FIDs/transients as
% (points x transients) and averaged spectra as (datasets x points)

rowsAreDatasets = false;

if isvector(d)
    d = d(:);
    return
end

if isempty(specDim)
    % Candidate lengths of the time/frequency dimension
    L = [];
    if isfield(MRS_struct, 'p')
        if isfield(MRS_struct.p, 'npoints') && numel(MRS_struct.p.npoints) >= ii
            L(end+1) = MRS_struct.p.npoints(ii);
        end
        if isfield(MRS_struct.p, 'ZeroFillTo') && numel(MRS_struct.p.ZeroFillTo) >= ii
            L(end+1) = MRS_struct.p.ZeroFillTo(ii);
        end
    end
    if isfield(MRS_struct, 'spec') && isfield(MRS_struct.spec, 'freq')
        L(end+1) = numel(MRS_struct.spec.freq);
    end
    in1 = ismember(size(d,1), L);
    in2 = ismember(size(d,2), L);
    if in1 && ~in2
        specDim = 1;
    elseif in2 && ~in1
        specDim = 2;
    elseif size(d,1) >= size(d,2)
        specDim = 1;
    else
        specDim = 2;
    end
end

if specDim == 2
    d = d.';
    rowsAreDatasets = true;
end

end


function domain = GuessDomain(d, label, isPath)
% Guess whether the data are FIDs (time domain) or spectra

if isPath && startsWith(label, 'fids')
    domain = 'time';
elseif isPath && startsWith(label, 'spec')
    domain = 'freq';
else
    % FIDs have their largest magnitude in the first few points
    [~, idx] = max(abs(d(:,1)));
    if idx <= 5
        domain = 'time';
    else
        domain = 'freq';
    end
end

end


function sw = GetSW(MRS_struct, ii, swOpt, isWater)
% Spectral width (Hz)

if ~isempty(swOpt)
    sw = swOpt;
elseif isWater && isfield(MRS_struct.p, 'sw_water') && numel(MRS_struct.p.sw_water) >= ii
    sw = MRS_struct.p.sw_water(ii);
elseif isfield(MRS_struct, 'p') && isfield(MRS_struct.p, 'sw') && numel(MRS_struct.p.sw) >= ii
    sw = MRS_struct.p.sw(ii);
else
    error('PlotMRSData:missingField', 'MRS_struct.p.sw(%d) does not exist. Supply ''SW'' or ''Freq''.', ii);
end

end


function nfft = GetZeroFillTo(MRS_struct, ii, sw)
% Zero-filled length (same rule as GannetLoad)

if isfield(MRS_struct.p, 'ZeroFillTo') && numel(MRS_struct.p.ZeroFillTo) >= ii
    nfft = MRS_struct.p.ZeroFillTo(ii);
else
    nfft = round(32768 / 2000 * sw);
end

end


function freq = GetFreqAxis(MRS_struct, ii, N, sw, opts)
% Use MRS_struct.spec.freq if it matches the data; otherwise create the ppm
% axis on the fly (same convention as GannetLoad)

if ~isempty(opts.Freq)
    if numel(opts.Freq) ~= N
        error('PlotMRSData:badFreq', 'The supplied ppm axis has %d points but the data have %d.', numel(opts.Freq), N);
    end
    freq = opts.Freq(:);
    return
end

if isempty(opts.F0) && isfield(MRS_struct, 'spec') && isfield(MRS_struct.spec, 'freq') ...
        && numel(MRS_struct.spec.freq) == N && isfield(MRS_struct.p, 'sw') ...
        && numel(MRS_struct.p.sw) >= ii && sw == MRS_struct.p.sw(ii)
    freq = MRS_struct.spec.freq(:);
    return
end

if ~isfield(MRS_struct.p, 'LarmorFreq') || numel(MRS_struct.p.LarmorFreq) < ii
    error('PlotMRSData:missingField', 'MRS_struct.p.LarmorFreq(%d) does not exist. Supply ''Freq''.', ii);
end

if ~isempty(opts.F0)
    F0 = opts.F0;
elseif isfield(MRS_struct.p, 'phantom') && MRS_struct.p.phantom
    F0 = 4.8;
else
    F0 = 4.68;
end

freqRange = sw / MRS_struct.p.LarmorFreq(ii);
freq = (N + 1 - (1:N)') / N * freqRange + F0 - freqRange/2;

end


function FormatAxes(ax, freqs, specs, xLim, yLim, isWater)
% Apply Gannet-style axis formatting and scale the y-axis. For non-water
% data, the y-axis is scaled to the signal between 0.5 and 4.1 ppm (upper
% limit 15% above the maximum) so that residual water does not dominate the
% plot; for water data, it is scaled to the data within xLim

if isempty(xLim)
    allFreq = cell2mat(cellfun(@(x) x(:), freqs, 'UniformOutput', false));
    xLim = [min(allFreq) max(allFreq)];
end
set(ax, 'XDir', 'reverse', 'XLim', sort(xLim), 'Box', 'off', 'TickDir', 'out');
xlabel(ax, 'ppm', 'FontSize', 16);

if isempty(yLim)
    if isWater
        yRefRange = xLim;
    else
        yRefRange = [0.5 4.1];
    end
    [yMin, yMax] = SignalRange(freqs, specs, yRefRange);
    if ~isWater && ~(isfinite(yMin) && isfinite(yMax))
        % No data between 0.5 and 4.1 ppm: fall back to the data within xLim
        [yMin, yMax] = SignalRange(freqs, specs, xLim);
    end
    if isfinite(yMin) && isfinite(yMax) && yMax > yMin
        pad = 0.05 * (yMax - yMin);
        if isWater
            yTop = yMax + pad;
        else
            yTop = yMax + 0.15 * abs(yMax);
            if yTop <= yMax % e.g., yMax = 0
                yTop = yMax + pad;
            end
        end
        set(ax, 'YLim', [yMin - pad, yTop]);
    end
else
    set(ax, 'YLim', yLim);
end

end


function [yMin, yMax] = SignalRange(freqs, specs, ppmRange)
% Minimum and maximum signal within a ppm range across all inputs

yMin = Inf;
yMax = -Inf;
for mm = 1:numel(specs)
    inRange = freqs{mm} >= min(ppmRange) & freqs{mm} <= max(ppmRange);
    if any(inRange)
        vals = specs{mm}(inRange,:);
        yMin = min(yMin, min(vals(:)));
        yMax = max(yMax, max(vals(:)));
    end
end

end
