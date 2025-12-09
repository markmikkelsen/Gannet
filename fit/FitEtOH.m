function [MRS_struct, freqBounds, plotBounds, residPlot] = FitEtOH(MRS_struct, freq, DIFF, vox, target, ii, jj, kk)

freqBounds = find(freq <= 1.8 & freq >= 0.6);
plotBounds = find(freq <= 1.9 & freq >= 0.4);

maxinEtOH = max(real(DIFF(ii,freqBounds)));
grad_points = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
LinearInit = grad_points ./ abs(freq(1) - freq(2));

LorentzModelInit = [maxinEtOH 1.11 1/500 ...
                    maxinEtOH 1.23 1/500 ...
                   -LinearInit 0];
LorentzModelInit([1 4 7]) = LorentzModelInit([1 4 7]) / maxinEtOH; % Scale initial conditions to avoid warnings about numerical underflow

lb = [0             1.11-0.01 1/700 0             1.23-0.01 1/700 -40*maxinEtOH -2e3*maxinEtOH];
ub = [100*maxinEtOH 1.11+0.01 1/300 100*maxinEtOH 1.23+0.01 1/300  40*maxinEtOH  1e3*maxinEtOH];

lb([1 4 7]) = lb([1 4 7]) / maxinEtOH;
ub([1 4 7]) = ub([1 4 7]) / maxinEtOH;

% Down-weight co-edited Lac signal by including observation
% weights in nonlinear regression
w = ones(size(DIFF(ii,freqBounds)));
residfreq = freq(freqBounds);
LacRange = residfreq >= 1.29 & residfreq <= 1.51;
weightRange = LacRange;
w(weightRange) = 0.001;

% Weighted least-squares model fitting
LorentzModelInit = lsqcurvefit(@EtOHModel, LorentzModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinEtOH, lb, ub, lsqopts);
modelFun_w = @(x,freq) sqrt(w) .* EtOHModel(x,freq); % add weights to the model
[LorentzModelParam, resid] = nlinfit(freq(freqBounds), sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinEtOH, ...
                                modelFun_w, LorentzModelInit, nlinopts); % add weights to the data
[~, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinEtOH, ...
                    @EtOHModel, LorentzModelParam, nlinopts); % re-run for residuals for output figure

% Rescale fit parameters and residuals
LorentzModelParam([1 4 7 8]) = LorentzModelParam([1 4 7 8]) * maxinEtOH;
resid = resid * maxinEtOH;
residPlot = residPlot * maxinEtOH;

EtOHheight = max(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqBounds)));
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / EtOHheight;
MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = sum(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqBounds))) * abs(freq(1) - freq(2));

MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii)         = (LorentzModelParam(3) + LorentzModelParam(6)) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = LorentzModelParam;
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:)      = resid;

% Calculate SNR of EtOH signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).EtOH.SNR(ii) = abs(EtOHheight) / noiseSigma_DIFF;
