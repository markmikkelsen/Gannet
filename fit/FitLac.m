function [MRS_struct, modelFit] = FitLac(MRS_struct, freq, spec, vox, ...
    target, ii, jj, kk, baseline, ~, ~, lsqnlinopts)

freqBounds = find(freq <= 1.8 & freq >= 0.5);
plotBounds = find(freq <= 2.12 & freq >= 0);

maxinLac = max(real(spec(freqBounds)));

                  % amplitude  % width  % freq
modelParamInit = [maxinLac     0.02     1.28 ...
                  maxinLac     0.02     1.33 ...
                  maxinLac/2  -50       1.18 ...
                  maxinLac/2  -500];
lb = [0      0  1.28-0.02 ...
      0      0  1.33-0.02 ...
      0  -2000  1.18-0.03 ...
      0  -2000];
ub = [2*maxinLac  0.2  1.28+0.02 ...
      2*maxinLac  0.2  1.33+0.02 ...
      maxinLac   -10   1.18+0.03 ...
      maxinLac   -10];

% Scale initial conditions and bounds to avoid warnings about numerical underflow
amplParams = [1 4 7 10];
modelParamInit(amplParams) = modelParamInit(amplParams) / maxinLac;
lb(amplParams) = lb(amplParams) / maxinLac;
ub(amplParams) = ub(amplParams) / maxinLac;

% Model fitting with fixed baseline
[modelParam, resid, h_tmp] = FitSignalModel(@LacModel_noBaseline, ... % model
                                freq(freqBounds), ... % freq
                                real(spec(freqBounds)) / maxinLac, ... % spec
                                baseline(freqBounds) / maxinLac, ... % baseline
                                modelParamInit, ... % beta0
                                lb, ...
                                ub, ...
                                lsqnlinopts);

% Rescale fit parameters and residuals
modelParam(amplParams) = modelParam(amplParams) * maxinLac;
resid = resid * maxinLac;

% Peak fits
[P1, P2, P3, P4] = deal(modelParam);
P1(amplParams(2:4))     = 0;
P2(amplParams([1 3:4])) = 0;
P3(amplParams([1:2 4])) = 0;
P4(amplParams(1:3))     = 0;

% Plot individual Gaussians
hold on;
plot(freq(freqBounds), LacModel_noBaseline(P1, freq(freqBounds)) / maxinLac);
plot(freq(freqBounds), LacModel_noBaseline(P2, freq(freqBounds)) / maxinLac);
plot(freq(freqBounds), LacModel_noBaseline(P3, freq(freqBounds)) / maxinLac);
plot(freq(freqBounds), LacModel_noBaseline(P4, freq(freqBounds)) / maxinLac);
legend({'data','model + baseline','baseline','residual', ...
    'Lac+ (1)','Lac+ (2)','BHB+','MM_{1.43}'}, ...
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
exportgraphics(h_tmp, fullfile(out_dir, [fname '_Lac_model_fit.png']), "Resolution", 300);
close(h_tmp);

MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = modelParam;
LacModelParam = modelParam;
LacModelParam(amplParams(3:4)) = 0;
MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = pi * (LacModelParam(amplParams(1)) * LacModelParam(2) + ...
                                                    LacModelParam(amplParams(2)) * LacModelParam(5)); % NB: this is Lac+
LacHeight = max(LacModel_noBaseline(modelParam, freq(freqBounds)));
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / LacHeight;
MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = 2 * (LacModelParam(2) + LacModelParam(5)) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;

% Calculate SNR of Lac+ signal
noiseSigma_DIFF = CalcNoise(freq, spec);
MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(LacHeight) / noiseSigma_DIFF;

modelFit.modelParam.full = modelParam;
modelFit.freqBounds      = freqBounds;
modelFit.plotBounds      = plotBounds;
modelFit.residPlot       = resid;

modelFit.modelParam.peak1 = P1;
modelFit.modelParam.peak2 = P2;
modelFit.modelParam.peak3 = P3;
modelFit.modelParam.peak4 = P4;
