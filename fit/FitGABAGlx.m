function [MRS_struct, modelFit] = FitGABAGlx(MRS_struct, freq, DIFF, vox, ...
    ~, ii, ~, kk, baseline, ~, ~, lsqnlinopts)

freqBounds = find(freq <= 4.1 & freq >= 2.79);
plotBounds = find(freq <= 4.2 & freq >= 2.7);

GABAbounds = freq <= 3.2 & freq >= 2.79;
Glxbounds  = freq <= 4.1 & freq >= 3.4;

maxinGABA   = max(real(DIFF(ii,GABAbounds)));
maxinGlx    = max(real(DIFF(ii,Glxbounds)));
grad_points = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
LinearInit  = grad_points ./ abs(freq(1) - freq(2));

GaussModelInit = [maxinGlx -700 3.71 maxinGlx -700 3.79 maxinGABA -90 3.02 -LinearInit 0 0];
GaussModelInit([1 4 7 10]) = GaussModelInit([1 4 7 10]) / maxinGlx; % Scale initial conditions to avoid warnings about numerical underflow

lb = [-4000*maxinGlx -1000 3.71-0.02 -4000*maxinGlx -1000 3.79-0.02 -4000*maxinGABA -200 3.02-0.05 -40*maxinGABA -2000*maxinGABA -2000*maxinGABA];
ub = [4000*maxinGlx  -40   3.71+0.02  4000*maxinGlx -40   3.79+0.02  4000*maxinGABA -40  3.02+0.05  40*maxinGABA  1000*maxinGABA  1000*maxinGABA];
lb([1 4 7 10]) = lb([1 4 7 10]) / maxinGlx;
ub([1 4 7 10]) = ub([1 4 7 10]) / maxinGlx;

% Down-weight Cho subtraction artifact and (if HERMES)
% signals downfield of Glx by including observation weights
% in nonlinear regression; improves accuracy of peak
% fittings (MM: 170701 - thanks to Alex Craven of
% University of Bergen for this idea)
w = ones(size(DIFF(ii,freqBounds)));
residfreq = freq(freqBounds);
ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
GlxDownfieldRange = residfreq >= 3.9 & residfreq <= 4.2;
if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor, {'Philips','Philips_data','Philips_raw'}))
    weightRange = ChoRange | GlxDownfieldRange;
else
    weightRange = ChoRange;
end
w(weightRange) = 0.001;

% Weighted least-squares model fitting
% GaussModelInit = lsqcurvefit(@GABAGlxModel, GaussModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGlx, lb, ub, lsqopts);
% modelFun_w = @(x,freq) sqrt(w) .* GABAGlxModel(x,freq); % add weights to the model
% [modelParam, resid] = nlinfit(freq(freqBounds), sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinGlx, ...
%                                 modelFun_w, GaussModelInit, nlinopts); % add weights to the data
% [~, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGlx, ...
%                     @GABAGlxModel, modelParam, nlinopts); % re-run for residuals for output figure

% Model fitting with predetermined baseline
GABAGlxModel_noBaseline_w = @(x,freq) sqrt(w).' .* GABAGlxModel_noBaseline(x,freq); % add weights to the model
% Model fitting with predetermined baseline
[modelParam, resid, h_tmp] = FitSignalModel(GABAGlxModel_noBaseline_w, ... % weighted model
                                freq(freqBounds), ... % freq
                                sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinGlx, ... % data
                                sqrt(w) .* baseline(freqBounds) / maxinGlx, ... % baseline
                                GaussModelInit(1:end-3), ... % beta0
                                lb(1:end-3), ...
                                ub(1:end-3), ...
                                lsqnlinopts);

% Rescale fit parameters and residuals
amplParams = [1 4 7];
modelParam(amplParams) = modelParam(amplParams) * maxinGlx;
resid     = resid * maxinGlx;
% residPlot = residPlot * maxinGlx;
residPlot = resid;

% Range to determine residuals for GABA and Glx
residGABA = resid(residfreq <= 3.55 & residfreq >= 2.79);
residGlx  = resid(residfreq <= 4.10 & residfreq >= 3.45);

% Gaussians
[Gauss1, Gauss2, Gauss3] = deal(modelParam);
Gauss1(amplParams([2 3])) = 0;
Gauss2(amplParams([1 3])) = 0;
Gauss3(amplParams([1 2])) = 0;

% Plot individual Gaussians
% figure(h_tmp);
hold on;
plot(freq(freqBounds), GABAGlxModel_noBaseline(Gauss1, freq(freqBounds)) / maxinGlx);
plot(freq(freqBounds), GABAGlxModel_noBaseline(Gauss2, freq(freqBounds)) / maxinGlx);
plot(freq(freqBounds), GABAGlxModel_noBaseline(Gauss3, freq(freqBounds)) / maxinGlx);
legend({'data','model + baseline','baseline','residual', ...
    'gauss1','gauss2','gauss3'}, ...
    'Box','off','Location','best');
hold off;
drawnow;

out_dir = fullfile(pwd, 'Gannet_model_output');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

[~,fname,ext] = fileparts(MRS_struct.metabfile{ii});
if strcmpi(ext, '.gz')
    fname(end-3:end) = [];
end
exportgraphics(h_tmp, fullfile(out_dir, [fname '_GABAGlx_model_fit.png']), "Resolution", 300);

% GABA fitting output
MRS_struct.out.(vox{kk}).GABA.Area(ii) = modelParam(7) ./ sqrt(-modelParam(8)) * sqrt(pi);
GABAheight = modelParam(7);
MRS_struct.out.(vox{kk}).GABA.FitError(ii) = 100 * std(residGABA) / GABAheight;
sigma = sqrt(1/(2*(abs(modelParam(8)))));
MRS_struct.out.(vox{kk}).GABA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).GABA.Resid(ii,:) = residGABA;

% Calculate SNR of GABA signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).GABA.SNR(ii) = abs(GABAheight) / noiseSigma_DIFF;

% MM (200728)
% MRS_struct.out.(vox{kk}).GABA.FitError2(ii) = sqrt(mean(residGABA.^2)) / (0.5*noiseSigma_DIFF);

% Glx fitting output
MRS_struct.out.(vox{kk}).Glx.Area(ii) = (modelParam(1) / sqrt(-modelParam(2)) * sqrt(pi)) + ...
                                            (modelParam(4) / sqrt(-modelParam(5)) * sqrt(pi));
Glxheight = max(modelParam([1,4]));
MRS_struct.out.(vox{kk}).Glx.FitError(ii) = 100 * std(residGlx) / Glxheight;
sigma = sqrt(1/(2*(abs(modelParam(2))))) + sqrt(1/(2*(abs(modelParam(5)))));
MRS_struct.out.(vox{kk}).Glx.FWHM(ii) = abs(2 * MRS_struct.p.LarmorFreq(ii) * sigma);
MRS_struct.out.(vox{kk}).Glx.ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).Glx.Resid(ii,:) = residGlx;

% Calculate SNR of Glx signal
MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight) / noiseSigma_DIFF;

% MM (200728)
% MRS_struct.out.(vox{kk}).Glx.FitError2(ii) = sqrt(mean(residGlx.^2)) / (0.5*noiseSigma_DIFF);

modelFit.modelParam        = modelParam;
modelFit.freqBounds        = freqBounds;
modelFit.plotBounds        = plotBounds;
modelFit.residPlot         = residPlot;
modelFit.weightRange       = weightRange;
modelFit.ChoRange          = ChoRange;
modelFit.GlxDownfieldRange = GlxDownfieldRange;

modelFit.ampl.Gauss1 = Gauss1;
modelFit.ampl.Gauss2 = Gauss2;
modelFit.ampl.Gauss3 = Gauss3;
