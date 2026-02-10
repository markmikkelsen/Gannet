function PaperPlot(MRS_struct, varargin)
% PaperPlot(MRS_struct, varargin)
%
% This function will plot the difference spectra saved in MRS_struct. The
% corresponding model fits can also be plotted. Users can choose to plot a
% single spectrum, a select number of spectra, or all spectra. Multiple
% spectra will be overlaid in the same figure. If data were acquired with
% HERMES, then each Hadamard-combined difference spectrum will be plotted
% in a separate subplot.
%
% If GannetCoRegister was run, users also have the option to plot an
% exemplary voxel mask co-registered to a corresponding structural image.
%
% To export plots at publication quality, consider using PaperPlot with
% Yair Altman's excellent export_fig toolbox
% (https://github.com/altmany/export_fig).
%
% Inputs:
%   MRS_struct: Structure output from GannetFit (required).
%   varargin:   Optional inputs (entered as parameter-value pairs).
%           target:     (For HERMES data only.) Choose a single target
%                       metabolite to plot, entered as a string. Default is
%                       plotting of difference spectra for all target
%                       metabolites.
%           specNum:    Spectra to plot, entered as a scalar or vector. All
%                       spectra are plotted by default.
%           freqLim:    Limits of ppm axis, entered as a two-element
%                       vector. Default is [0.5 4.5].
%           signalLim:  Limits of signal axis, entered as a two-element
%                       vector. Default is an empty vector (automatic
%                       scaling; recommended).
%           plotModel:  Plot signal model fit(s), entered as a logical.
%                       Default is false.
%           plotResid:  If plotModel is true, also plot the model residuals,
%                       entered as a logical. Default is true.
%           plotAvg:    Plot the group-average spectrum, entered as a
%                       logical. Default is false.
%           plotStd:    If plotAvg is true, also show the +/- 1 standard
%                       deviation, entered as a logical. Default is false.
%           plotCI:     If plotAvg is true, also show the 95% confidence
%                       interval, entered as a logical. Default is false.
%           plotVoxMask:    Plot an exemplary voxel mask co-registered to
%                           the respective structural image
%                           (GannetCoRegister.m must have been run).
%           voxNum:     Voxel mask to plot, entered as a scalar. Default is
%                       the voxel from the first dataset in MRS_struct.
%           plotPreAlign:   Plot the difference spectra before frequency
%                           and phase alignment, entered as a logical.
%                           Default is false.
%
% Examples:
%   PaperPlot(MRS_struct, 'specNum', [1 3 4]);
%       This will plot the 1st, 3rd and 4th difference spectra in
%       MRS_struct along with the model fits of the peak(s) specified in
%       MRS_struct.p.target.
%
%   PaperPlot(MRS_struct, 'freqLim', [2.5 3.5], 'plotModel', true);
%       This will plot all difference spectra in MRS_struct with the model
%       fits of the peak(s) and limit the ppm axis from 2.5 to 3.5 ppm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Parse inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs','Not enough input arguments.');
end

if ~isstruct(MRS_struct)
    fprintf('\n');
    error('The first input argument must be a structure, but received %s.', class(MRS_struct));
end

% Set some defaults
vox = MRS_struct.p.vox;
if ~MRS_struct.p.PRIAM
    vox = vox(1);
end
defaultTarget       = MRS_struct.p.target;
defaultnSpec        = 1:length(MRS_struct.metabfile);
defaultFreqLim      = [0.5 4.5];
defaultSignalLim    = [];
defaultPlotModel    = false;
defaultPlotResid    = true;
defaultPlotAvg      = false;
defaultPlotStd      = false;
defaultPlotCI       = false;
expectedTargets     = {'GABAGlx', 'GSH', 'Lac', 'EtOH', 'GABA', 'Glx'};
defaultPlotVoxMask  = false;
defaultnVox         = 1;
defaultPlotPreAlign = false;
grey                = [0.6 0.6 0.6];
shading             = 0.3;

% Parse input arguments
p = inputParser;
p.CaseSensitive = false;
p.addParameter('target', defaultTarget, @(x) any(validatestring(x,expectedTargets)));
p.addParameter('specNum', defaultnSpec, @(x) isvector(x));
p.addParameter('freqLim', defaultFreqLim, @(x) isvector(x) && numel(x) == 2);
p.addParameter('signalLim', defaultSignalLim, @(x) isvector(x) && numel(x) == 2);
p.addParameter('plotModel', defaultPlotModel, @(x) islogical(x));
p.addParameter('plotResid', defaultPlotResid, @(x) islogical(x));
p.addParameter('plotAvg', defaultPlotAvg, @(x) islogical(x));
p.addParameter('plotStd', defaultPlotStd, @(x) islogical(x));
p.addParameter('plotCI', defaultPlotCI, @(x) islogical(x));
p.addParameter('plotVoxMask', defaultPlotVoxMask, @(x) islogical(x));
p.addParameter('voxNum', defaultnVox, @(x) isscalar(x));
p.addParameter('plotPreAlign', defaultPlotPreAlign, @(x) islogical(x));
p.parse(varargin{:});

target = p.Results.target;
if ischar(target)
    target = {target};
end
specNum      = p.Results.specNum;
freqLim      = p.Results.freqLim;
signalLim    = p.Results.signalLim;
plotModel    = p.Results.plotModel;
plotResid    = p.Results.plotResid;
plotAvg      = p.Results.plotAvg;
plotStd      = p.Results.plotStd;
plotCI       = p.Results.plotCI;
plotVoxMask  = p.Results.plotVoxMask;
voxNum       = p.Results.voxNum;
plotPreAlign = p.Results.plotPreAlign;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq = MRS_struct.spec.freq;

if isfield(MRS_struct.spec.(vox{1}).(target{1}), 'diff_scaled')
    diff = 'diff_scaled';
    fprintf('\nNB: Spectra are normalized to the amplitude of the respective modeled unsuppressed water reference signal.\n\n');
else
    diff = 'diff';
end

if plotPreAlign
    diff = 'diff_noalign';
    fprintf('NB: Plotting the pre-aligned difference spectra.\n\n');
end

if plotVoxMask && ~isfield(MRS_struct, 'mask')
    error('GannetCoRegister.m has not been run. Cannot plot voxel mask image.');
end

for ii = 1:length(vox)

    H = figure(199+ii);
    scr_sz = get(0, 'ScreenSize');
    fig_w = 1000;
    if length(target) > 1
        fig_h = 1000;
    else
        fig_h = 500;
    end
    set(H, 'Color', 'w', 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
    set(H, 'Name', 'PaperPlot (Spectra)', 'Tag', 'PaperPlot', 'NumberTitle', 'off');
    clf;

    if isfield(MRS_struct.out.(vox{ii}), 'water')
        scaleFactor = MRS_struct.out.(vox{ii}).water.ModelParam(specNum,1);
    else
        scaleFactor = ones(1,length(specNum));
    end

    for jj = 1:length(target)

        if length(target) > 1
            H = subplot(length(target),1,jj);
        else
            H = gca;
        end

        switch target{jj}
            case 'GABA'
                if length(target) == 3 && all(ismember(target, {'EtOH', 'GABA', 'GSH'}))
                    modelFreq = freq(freq <= 3.55 & freq >= 2.6);
                    residInd  = freq <= 3.55 & freq >= 2.6;
                else
                    modelFreq = freq(freq <= 3.55 & freq >= 2.79);
                    residInd  = freq <= 3.55 & freq >= 2.79;
                end
                model         = @GaussModel;
                baselineFreq  = freq <= 3.5 & freq >= 3.4;
            case 'GABAGlx'
                modelFreq     = freq(freq <= 4.1 & freq >= 2.79);
                model         = @GABAGlxModel_noBaseline;
                residInd      = freq <= 4.1 & freq >= 2.79;
                baselineFreq  = freq <= 3.5 & freq >= 3.4;
            case 'GSH'
                modelFreq     = freq(freq <= 3.5 & freq >= 2.25);
                model         = @EightGaussModel_noBaseline;
                residInd      = freq <= 3.5 & freq >= 2.25;
                baselineFreq  = freq <= 1.8 & freq >= 1.7;
            case 'Lac'
                modelFreq     = freq(freq <= 1.8 & freq >= 0.5);
                model         = @LacModel;
                residInd      = freq <= 1.8 & freq >= 0.5;
                baselineFreq  = freq <= 0.25 & freq >= -0.25;
            case 'EtOH'
                modelFreq     = freq(freq <= 1.8 & freq >= 0.6);
                model         = @EtOHModel;
                residInd      = freq <= 1.8 & freq >= 0.6;
                baselineFreq  = freq <= 0.25 & freq >= -0.25;
        end

        % Demean baseline
        baseMean = repmat(mean(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum, baselineFreq)),2), ...
            [1 size(MRS_struct.spec.(vox{ii}).(target{jj}).(diff),2)]);
        baseMeanResid = repmat(mean(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum, baselineFreq)),2), [1 length(modelFreq)]);

        if numel(specNum) > 1 && plotAvg

            % Find mean, std and 95% CI
            mu       = mean(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum,:)) - baseMean,1);
            sigma    = std(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum,:)) - baseMean,[],1);
            stdErr   = sigma / sqrt(numel(specNum));
            UB.sigma = mu + sigma;
            LB.sigma = mu - sigma;
            UB.ci    = mu + 1.96 * stdErr;
            LB.ci    = mu - 1.96 * stdErr;

            hold on;
            if plotStd && plotCI
                patch([freq fliplr(freq)], [UB.sigma fliplr(LB.sigma)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
                patch([freq fliplr(freq)], [UB.ci fliplr(LB.ci)], 1, 'FaceColor', (grey-0.4)+(1-(grey-0.4))*(1-shading), 'EdgeColor', 'none');
            elseif plotStd
                patch([freq fliplr(freq)], [UB.sigma fliplr(LB.sigma)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            elseif plotCI
                patch([freq fliplr(freq)], [UB.ci fliplr(LB.ci)], 1, 'FaceColor', grey+(1-grey)*(1-shading), 'EdgeColor', 'none');
            end
            h = plot(freq, mu, 'Color', [0 0 0], 'LineWidth', 1);
            hold off;

        else

            hold on;
            for kk = 1:numel(specNum)
                if plotModel
                    if strcmp(target{jj}, 'GABAGlx')
                        h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), ...
                            modelFreq, model(MRS_struct.out.(vox{ii}).GABA.ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk) - baseMean(kk,1), ...
                            'LineWidth', 1);
                    else
                        h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), ...
                            modelFreq, model(MRS_struct.out.(vox{ii}).(target{jj}).ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk) - baseMean(kk,1), 'LineWidth', 1);
                    end
                    h(1,kk).Color = [0 0 0];
                    h(2,kk).Color = [1 0 0];
                else
                    h(:,kk) = plot(freq, real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),:)) - baseMean(kk,:), 'Color', [0 0 0], 'LineWidth', 1);
                end
                if plotResid && plotModel
                    if strcmp(target{jj}, 'GABAGlx')
                        resid = (real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),residInd)) - ...
                            (model(MRS_struct.out.(vox{ii}).GABA.ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk))) - baseMeanResid(kk,:);
                    else
                        resid = (real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk),residInd)) - ...
                            (model(MRS_struct.out.(vox{ii}).(target{jj}).ModelParam(specNum(kk),:),modelFreq) ./ scaleFactor(kk))) - baseMeanResid(kk,:);
                    end
                    dataMin = min(real(MRS_struct.spec.(vox{ii}).(target{jj}).(diff)(specNum(kk), residInd)),[],2);
                    resid = resid + dataMin - 1.5*max(resid);
                    plot(modelFreq, resid, 'Color', [0 0 0], 'LineWidth', 1);
                end
            end
            hold off;

        end

        % Set YLim
        if isempty(signalLim)
            switch target{jj}
                case {'GABAGlx', 'GABA', 'Glx'}
                    if MRS_struct.p.phantom
                        peakRange = freq <= 4.25 & freq >= 1.0;
                    else
                        peakRange = freq <= 4.1 & freq >= 2.26;
                    end
                case 'GSH'
                    peakRange = freq <= 3.5 & freq >= 0.5;
                case {'Lac', 'EtOH'}
                    peakRange = freq <= 3 & freq >= 0.5;
            end
            if plotStd && plotCI
                yAxisMax = max([UB.sigma(peakRange) UB.ci(peakRange)]);
                yAxisMin = min([LB.sigma(peakRange) LB.ci(peakRange)]);
            elseif plotStd
                yAxisMax = max(UB.sigma(peakRange));
                yAxisMin = min(LB.sigma(peakRange));
            elseif plotCI
                yAxisMax = max(UB.ci(peakRange));
                yAxisMin = min(LB.ci(peakRange));
            else
                for kk = 1:size(h,2)
                    maxPeakHeight(kk) = max(h(1,kk).YData(peakRange)); %#ok<*AGROW>
                    minPeakHeight(kk) = min(h(1,kk).YData(peakRange));
                end
                yAxisMax = max(maxPeakHeight);
                yAxisMin = min(minPeakHeight);
            end
            yRange = abs(yAxisMax - yAxisMin);
            yAxisMax = yAxisMax + 0.1*yRange;
            if any(strcmp(target{jj}, {'GABAGlx', 'GABA', 'Glx'}))
                if MRS_struct.p.phantom
                    yAxisMin = yAxisMin - 0.15*abs(yAxisMin);
                else
                    yAxisMin = yAxisMin - 0.3*yRange;
                end
            else
                yAxisMin = yAxisMin - 0.1*yRange;
            end
            signalLim = [yAxisMin yAxisMax];
            set(H, 'YLim', signalLim);
            signalLim = [];
        else
            set(H, 'YLim', signalLim(jj,:));
        end

        set(gca, 'TickDir', 'out', 'XLim', freqLim, 'XDir', 'reverse', 'Box', 'off', ...
            'FontSize', 20, 'LineWidth', 1, 'XColor', [0 0 0], 'YColor', [0 0 0]);
        set(get(gca,'YAxis'),'Visible','off');
        xlabel('ppm', 'FontWeight', 'bold', 'FontSize', 28, 'Color', [0 0 0]);

    end

    if plotVoxMask

        H2 = figure(200+ii);
        scr_sz = get(0, 'ScreenSize');
        fig_w = 1000;
        fig_h = round(fig_w / (size(MRS_struct.mask.vox1.img{voxNum},2) / size(MRS_struct.mask.vox1.img{voxNum},1)));
        set(H2, 'Color', 'w', 'Position', [(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(H2, 'Name', 'PaperPlot (Voxel Mask)', 'Tag', 'PaperPlot', 'NumberTitle', 'off');
        clf;

        axes('Position', [0 0 1 1]);
        imagesc(MRS_struct.mask.(vox{ii}).img{voxNum});
        colormap('gray');
        img = MRS_struct.mask.(vox{ii}).img{voxNum}(:);
        caxis([0 mean(img(img > 0.01)) + 3*std(img(img > 0.01))]); %#ok<CAXIS>
        axis equal;
        axis tight;
        axis off;
        text(10, size(MRS_struct.mask.(vox{ii}).img{voxNum},1)/2, 'L', 'Color', [1 1 1], 'FontSize', 20);
        text(size(MRS_struct.mask.(vox{ii}).img{voxNum},2) - 20, size(MRS_struct.mask.(vox{ii}).img{voxNum},1)/2, 'R', 'Color', [1 1 1], 'FontSize', 20);

        set(findall(H2,'-property','FontName'),'FontName','Arial');

    end

    set(findall(H,'-property','FontName'),'FontName','Arial');

end
