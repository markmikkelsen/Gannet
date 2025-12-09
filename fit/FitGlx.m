function [MRS_struct, freqBounds, plotBounds, residPlot] = FitGlx(MRS_struct, freq, DIFF, vox, target, ii, jj, kk, lsqopts, nlinopts)

freqBounds = find(freq <= 4.1 & freq >= 3.45);
plotBounds = find(freq <= 4.5 & freq >= 3);

maxinGlx    = max(real(DIFF(ii,freqBounds)));
grad_points = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
LinearInit  = grad_points ./ abs(freq(1) - freq(2));
constInit   = (real(DIFF(ii,freqBounds(end))) + real(DIFF(ii,freqBounds(1))))./2;

GaussModelInit = [maxinGlx -90 3.72 maxinGlx -90 3.77 -LinearInit constInit];
lb = [0 -200 3.72-0.01 0 -200 3.77-0.01 -40*maxinGlx -2000*maxinGlx];
ub = [4000*maxinGlx -40 3.72+0.01 4000*maxinGlx -40 3.77+0.01 40*maxinGlx 1000*maxinGlx];

GaussModelInit = lsqcurvefit(@DoubleGaussModel, GaussModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)), lb, ub, lsqopts);
[GaussModelParam, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)), @DoubleGaussModel, GaussModelInit, nlinopts);

Glxheight = max(GaussModelParam([1,4]));
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(residPlot) / Glxheight;
MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = (GaussModelParam(1) / sqrt(-GaussModelParam(2)) * sqrt(pi)) + ...
                                                    (GaussModelParam(4) / sqrt(-GaussModelParam(5)) * sqrt(pi));
sigma = (sqrt(1/(2*(abs(GaussModelParam(2)))))) + (sqrt(1/(2*(abs(GaussModelParam(5))))));
MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * sigma);
MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = residPlot;

% Calculate SNR of Glx signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight) / noiseSigma_DIFF;
