function [MRS_struct, modelFit] = FitEtOH(MRS_struct, freq, spec, vox, ...
    target, ii, jj, kk, baseline, ~, ~, lsqnlinopts)

freqBounds = find(freq <= 1.8 & freq >= 0.6);
plotBounds = find(freq <= 1.9 & freq >= 0.4);

maxinEtOH = max(real(spec(freqBounds)));
grad_points = (real(spec(freqBounds(end))) - real(spec(freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
LinearInit = grad_points ./ abs(freq(1) - freq(2));

LorentzModelInit = [maxinEtOH  1.11  1/500 ...
                    maxinEtOH  1.23  1/500 ...
                   -LinearInit  0];
LorentzModelInit([1 4 7]) = LorentzModelInit([1 4 7]) / maxinEtOH; % Scale initial conditions to avoid warnings about numerical underflow

lb = [0  1.11-0.01  1/700 ...
      0  1.23-0.01  1/700 ...
      -40*maxinEtOH  -2e3*maxinEtOH];
ub = [100*maxinEtOH  1.11+0.01  1/300 ...
      100*maxinEtOH  1.23+0.01  1/300 ...
      40*maxinEtOH  1e3*maxinEtOH];

lb([1 4 7]) = lb([1 4 7]) / maxinEtOH;
ub([1 4 7]) = ub([1 4 7]) / maxinEtOH;

% Down-weight co-edited Lac signal by including observation
% weights in nonlinear regression
w = ones(size(spec(freqBounds)));
residfreq = freq(freqBounds);
LacRange = residfreq >= 1.29 & residfreq <= 1.51;
weightRange = LacRange;
w(weightRange) = 0.001;

% Weighted least-squares model fitting with fixed baseline
EtOHModel_noBaseline_w = @(x,freq) sqrt(w).' .* EtOHModel_noBaseline(x,freq); % add weights to the model
[modelParam, resid, h_tmp] = FitSignalModel(EtOHModel_noBaseline_w, ... % weighted model
                                freq(freqBounds), ... % freq
                                sqrt(w) .* real(spec(freqBounds)) / maxinGlx, ... % spec
                                sqrt(w) .* baseline(freqBounds) / maxinGlx, ... % baseline
                                modelParamInit, ... % beta0
                                lb, ...
                                ub, ...
                                lsqnlinopts);
% Re-run for residuals for output figure
[~, residPlot] = FitSignalModel(@EtOHModel_noBaseline, ... % weighted model (@ is needed here to avoid an error)
                    freq(freqBounds), ... % freq
                    real(spec(freqBounds)) / maxinGlx, ... % data
                    baseline(freqBounds) / maxinGlx, ... % baseline
                    modelParam, ... % beta0
                    lb, ...
                    ub, ...
                    lsqnlinopts);

% Rescale fit parameters and residuals
modelParam([1 4 7 8]) = modelParam([1 4 7 8]) * maxinEtOH;
resid = resid * maxinEtOH;
residPlot = residPlot * maxinEtOH;

EtOHheight = max(EtOHModel([modelParam(1:end-2) 0 0],freq(freqBounds)));
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / EtOHheight;
MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = sum(EtOHModel([modelParam(1:end-2) 0 0],freq(freqBounds))) * abs(freq(1) - freq(2));

MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii)         = (modelParam(3) + modelParam(6)) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:)      = resid;

% Calculate SNR of EtOH signal
noiseSigma_DIFF = CalcNoise(freq, spec);
MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(EtOHheight) / noiseSigma_DIFF;

modelFit.modelParam.full = modelParam;
modelFit.freqBounds      = freqBounds;
modelFit.plotBounds      = plotBounds;
modelFit.residPlot       = residPlot;

modelFit.modelParam.peak1 = Lorentz1;
modelFit.modelParam.peak2 = Lorentz1;
