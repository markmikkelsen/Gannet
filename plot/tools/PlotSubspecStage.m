function MRS_struct = PlotSubspecStage(MRS_struct, data, ii, stage, stageTag, figNum, domain, reject)
% Plot the averaged ON/OFF (MEGA) or A–D (HERMES) subspectra at a given
% point in the pipeline using PlotMRSData. Transients are assigned to
% subspectra in the same way as in SignalAveraging. Called from GannetLoad
% when MRS_struct.p.plot_subspec = 1. The plot is always
% printed as a PNG to the subject-specific folder in Gannet_spec_plots
% (see GetSpecPlotDir), whether or not MRS_struct.p.hide = 1. The file
% name starts with a number from the plot counter MRS_struct.plot.num
% (reset for each dataset and voxel in GannetLoad and incremented here), so
% that the PNGs sort in pipeline order whichever processing steps are run.
% The name of the current voxel is taken from MRS_struct.plot.vox (also set
% in GannetLoad)
%
% Inputs:
%   data     - FIDs (points x transients) or spectra (points x transients)
%   ii       - dataset index
%   stage    - processing stage (used in the figure title)
%   stageTag - short stage label (used in the PNG file name)
%   figNum   - figure number to (re)use
%   domain   - 'time' (FIDs) or 'freq' (spectra)
%   reject   - logical vector of transients to exclude (optional)
%
% Output:
%   MRS_struct - with the plot counter (MRS_struct.plot.num) incremented
%
% Author: Mark Mikkelsen, Ph.D.
% Created using Claude Opus 5.5 (2026-11-04); human-inspected by the author

if MRS_struct.p.HERMES
    n_subspec  = 4;
    experiment = {'A','B','C','D'};
    subspecTag = 'A-D';
else
    n_subspec  = 2;
    experiment = {'ON','OFF'};
    subspecTag = 'ON_OFF';
    if ~MRS_struct.fids.ON_OFF(1) % OFF is first
        experiment = flip(experiment);
    end
end

if nargin < 8 || isempty(reject)
    reject = false(1, size(data,2));
end

% Average the transients of each subspectrum
subspec = cell(1, n_subspec);
for ll = 1:n_subspec
    if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
        ind = MRS_struct.fids.ON_OFF == abs(ll-2);
    else
        ind = ismember(1:size(data,2), ll:n_subspec:size(data,2));
    end
    subspec{ll} = mean(data(:, ind & ~reject), 2);
end

% Subject-specific plot folder and data file name (for the title and PNG)
[plotDir, fname] = GetSpecPlotDir(MRS_struct, ii);

% Apply the same line-broadening to FIDs as is applied to AllFramesFT so
% that all stages are comparable
if strcmp(domain, 'time')
    LB = MRS_struct.p.LB;
else
    LB = 0;
end

% Increment the plot counter (starts at 0 if it does not exist, e.g., when
% used outside GannetLoad)
if ~isfield(MRS_struct, 'plot') || ~isfield(MRS_struct.plot, 'num')
    MRS_struct.plot.num = 0;
end
MRS_struct.plot.num = MRS_struct.plot.num + 1;
plotNum = MRS_struct.plot.num;

% Current voxel (the first voxel if not set, e.g., when used outside
% GannetLoad)
if isfield(MRS_struct.plot, 'vox')
    voxName = MRS_struct.plot.vox;
else
    voxName = MRS_struct.p.vox{1};
end

% PNG file name (add the voxel name if there is more than one, i.e., PRIAM)
if MRS_struct.p.PRIAM && numel(MRS_struct.p.vox) > 1
    pngName = sprintf('%02d_%s_%s_%s_subspec_%s.png', plotNum, fname, voxName, subspecTag, stageTag);
else
    pngName = sprintf('%02d_%s_%s_subspec_%s.png', plotNum, fname, subspecTag, stageTag);
end

PlotMRSData(MRS_struct, subspec, 'ii', ii, 'Domain', domain, 'LB', LB, ...
    'Legend', experiment, 'FigNum', figNum, 'Title', {fname, [voxName ': ' stage]}, ...
    'OutDir', plotDir, 'FileName', pngName);
