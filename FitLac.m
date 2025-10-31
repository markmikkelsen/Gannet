function [MRS_struct, freqBounds, plotBounds, residPlot] = FitLac(MRS_struct, freq, DIFF, vox, ~, ii, ~, kk, lsqopts, nlinopts)

freqBounds = find(freq <= 1.8 & freq >= 0.5);
plotBounds = find(freq <= 2.12 & freq >= 0);

offset   = (mean(real(DIFF(ii,freqBounds(1:10))),2) + mean(real(DIFF(ii,freqBounds((end-9):end))),2)) / 2;
slope    = (mean(real(DIFF(ii,freqBounds(1:10))),2) - mean(real(DIFF(ii,freqBounds((end-9):end))),2)) / abs(freq(freqBounds(1)) - freq(freqBounds(end)));
maxinLac = max(real(DIFF(ii,freqBounds)));

LacModelInit = [maxinLac/2 -1000 1.31 ...
                maxinLac/4 -90   1.21 ...
                offset slope 0 ...
                maxinLac/2];
lb = [0 -4000 1.31-0.02 ...
      0 -100  1.24-0.02 ...
     -1 -1 -1 ...
      0];
ub = [maxinLac -500 1.31+0.02 ...
      maxinLac  0   1.24+0.02 ...
      1 1 1 ...
      maxinLac];

LacModelInit = lsqcurvefit(@LacModel, LacModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)), lb, ub,lsqopts);
[LacPlusModelParam, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)), @LacModel, LacModelInit, nlinopts);

MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:) = LacPlusModelParam;
%MMmodelParam         = LacPlusModelParam;
%MMmodelParam([4 7])  = 0;
LacModelParam         = LacPlusModelParam;
LacModelParam([4 10]) = 0;
BHBmodelParam         = LacPlusModelParam;
BHBmodelParam(1)      = 0;

MRS_struct.out.(vox{kk}).Lac.Area(ii) = sum(LacModel([LacPlusModelParam(1:3) zeros(1,7)], freq(freqBounds))) * abs(freq(1) - freq(2)); % NB: this is Lac+
modelHeight = max(LacModel([LacPlusModelParam(1:6) zeros(1,3) LacPlusModelParam(10)], freq(freqBounds)));
MRS_struct.out.(vox{kk}).Lac.FitError(ii) = 100 * std(residPlot) / modelHeight;
MRS_struct.out.(vox{kk}).Lac.FWHM(ii) = NaN; % MM (170818): Still need to calculate FWHM
MRS_struct.out.(vox{kk}).Lac.Resid(ii,:) = residPlot;

% Calculate SNR of Lac signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).Lac.SNR(ii) = abs(modelHeight) / noiseSigma_DIFF;
