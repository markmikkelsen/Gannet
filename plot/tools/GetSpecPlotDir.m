function [plotDir, fname] = GetSpecPlotDir(MRS_struct, ii, makeDir)
% Return the subject-specific folder for spectrum plots (PNGs) inside the
% Gannet_spec_plots folder, and the name of the data file it is named after:
%   non-BIDS: <current folder>/Gannet_spec_plots/<data file name>/
%   BIDS:     <BIDS root>/derivatives/Gannet_spec_plots/<data file name>/
% Both folders are created if they do not exist. If no data file name is
% available, the Gannet_spec_plots folder itself is returned
%
% Inputs:
%   MRS_struct - Gannet data structure
%   ii         - dataset index
%   makeDir    - create the folder if it does not exist (default: true)
%
% Outputs:
%   plotDir - full path of the folder
%   fname   - data file name without path and extension ('' if unknown)
%
% Author: Mark Mikkelsen, Ph.D.
% Created using Claude Opus 5.5 (2026-11-04); human-inspected by the author

if nargin < 3
    makeDir = true;
end

% Data file name (same convention as SavePDF)
fname = '';
if isfield(MRS_struct, 'metabfile') && iscell(MRS_struct.metabfile)
    isRDA = isfield(MRS_struct, 'p') && isfield(MRS_struct.p, 'vendor') ...
            && strcmp(MRS_struct.p.vendor, 'Siemens_rda');
    if isRDA && size(MRS_struct.metabfile,2) >= ii*2-1
        [~, fname] = fileparts(MRS_struct.metabfile{1,ii*2-1});
    elseif ~isRDA && size(MRS_struct.metabfile,2) >= ii
        [~, fname] = fileparts(MRS_struct.metabfile{1,ii});
    end
end

% Top-level folder for spectrum plots
if isfield(MRS_struct, 'p') && isfield(MRS_struct.p, 'bids') && MRS_struct.p.bids
    specPlotDir = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_spec_plots');
else
    specPlotDir = fullfile(pwd, 'Gannet_spec_plots');
end

plotDir = fullfile(specPlotDir, fname);
if makeDir && ~isfolder(plotDir)
    mkdir(plotDir);
end

end
