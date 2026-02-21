function [MRS_struct, modelFit] = FitGSH(MRS_struct, freq, DIFF, vox, ...
    target, ii, jj, kk, baseline, ~, ~, lsqnlinopts)

freqBounds = find(freq <= 3.5 & freq >= 2.1);
plotBounds = find(freq <= 4.2 & freq >= 1.75);

GSHbounds = freq <= 3.1 & freq >= 2.8;
AspBounds = freq <= 2.8 & freq >= 2.3;

maxinGSH = max(abs(real(DIFF(ii,GSHbounds))));
[maxinAsp, maxInd] = max(abs(real(DIFF(ii,AspBounds))));

DIFF_Asp = DIFF(ii,AspBounds);
s = sign(real(DIFF_Asp(maxInd)));
maxinAsp = s * maxinAsp;

if MRS_struct.p.HERMES
    s = -1;
else
    s = 1;
end

GSHgaussModel = @EightGaussModel_noBaseline;

if MRS_struct.p.TE(ii) < 100

                       % amplitude      % width  % freq  
    modelParamInit = [ maxinGSH         -350     2.95 ...
                       0                -300          ...
                       maxinGSH         -300     2.71 ...
                       maxinAsp*0.5     -1000         ...
                       maxinAsp         -500     2.56 ...
                       maxinAsp*0.5     -1000         ...
                       s*maxinAsp*0.15  -300     2.46 ...
                      -maxinGSH         -700     2.36];

    lb = [ 0               -2000  2.95-0.02 ...
          -2*maxinGSH      -5000            ...
           0               -5000  2.71-0.02 ...
           2*maxinAsp*0.5  -5000            ...
           2*maxinAsp      -5000  2.56-0.01 ...
           2*maxinAsp*0.5  -5000            ...
           0               -2000  2.46-0.01 ...
          -2*maxinGSH      -2000  2.36-0.02];

    ub = [ 2*maxinGSH     -125  2.95+0.02 ...
           2*maxinGSH     -300            ...
           1.5*maxinGSH   -300  2.71+0.02 ...
           0              -125            ...
           0              -125  2.56+0.01 ...
           0              -125            ...
          -1.25*maxinAsp  -300  2.46+0.01 ...
           0              -125  2.36+0.02];


else

                       % amplitude      % width  % freq  
    modelParamInit = [ maxinGSH         -350     2.95 ...
                      -maxinGSH         -300          ...
                       maxinGSH         -500     2.70 ...
                       maxinAsp*0.5     -1000         ...
                       maxinAsp         -1000    2.56 ...
                       maxinAsp*0.5     -1000         ...
                       s*maxinAsp*0.15  -1000    2.46 ...
                       maxinGSH         -1000    2.36];

    lb = [ 0           -2000  2.95-0.02 ...
          -2*maxinGSH  -5000            ...
           0           -5000  2.70-0.02 ...
           0           -5000            ...
           0           -5000  2.56-0.01 ...
           0           -5000            ...
           0           -2000  2.46-0.01 ...
           0           -2000  2.36-0.02];

    ub = [2*maxinGSH      -125  2.95+0.02 ...
          2*maxinGSH      -300            ...
          1.5*maxinGSH    -300  2.70+0.02 ...
          2*maxinAsp*0.5  -125            ...
          2*maxinAsp      -125  2.56+0.01 ...
          2*maxinAsp*0.5  -125            ...
          1.25*maxinAsp   -300  2.46+0.01 ...
          2*maxinGSH      -125  2.36+0.02];

end

% Scale initial conditions and bounds to avoid warnings about numerical underflow
amplParams = [1 4 6 9 11 14 16 19];
modelParamInit(amplParams) = modelParamInit(amplParams) / maxinGSH;
lb(amplParams) = lb(amplParams) / maxinGSH;
ub(amplParams) = ub(amplParams) / maxinGSH;

% Model fitting with fixed baseline
[modelParam, resid, h_tmp] = FitSignalModel(GSHgaussModel, ... % model
                                freq(freqBounds), ... % freq
                                real(DIFF(ii,freqBounds)) / maxinGSH, ... % data
                                baseline(freqBounds) / maxinGSH, ... % baseline
                                modelParamInit, ... % beta0
                                lb, ...
                                ub, ...
                                lsqnlinopts);

% Rescale fit parameters and residuals
modelParam(amplParams) = modelParam(amplParams) * maxinGSH;
resid = resid * maxinGSH;
% residPlot = residPlot * maxinGSH;
residPlot = resid;

% Range to determine residuals for GSH
residFreq = freq(freqBounds);
residFreqRange = residFreq <= 3.1 & residFreq >= 2.82;
residGSH = resid(residFreqRange);

% Gaussians
[Gauss1, Gauss2, Gauss3, Gauss4, ...
    Gauss5, Gauss6, Gauss7, Gauss8] = deal(modelParam);
Gauss1(amplParams(2:end))       = 0;
Gauss2(amplParams([1 3:end]))   = 0;
Gauss3(amplParams([1:2 4:end])) = 0;
Gauss4(amplParams([1:3 5:end])) = 0;
Gauss5(amplParams([1:4 6:end])) = 0;
Gauss6(amplParams([1:5 7:end])) = 0;
Gauss7(amplParams([1:6 8]))     = 0;
Gauss8(amplParams(1:end-1))     = 0;

% Plot individual Gaussians
hold on;
plot(freq(freqBounds), GSHgaussModel(Gauss1, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss2, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss3, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss4, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss5, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss6, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss7, freq(freqBounds)) / maxinGSH);
plot(freq(freqBounds), GSHgaussModel(Gauss8, freq(freqBounds)) / maxinGSH);
legend({'data','model + baseline','baseline','residual', ...
    'gauss1','gauss2','gauss3','gauss4','gauss5','gauss6','gauss7','gauss8'}, ...
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
% exportgraphics(h_tmp, fullfile(out_dir, [fname '_GSH_model_fit.png']), "Resolution", 300);
close(h_tmp);

MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = modelParam(1) ./ sqrt(-modelParam(2)) * sqrt(pi);
GSHheight = modelParam(1);
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(residGSH) / GSHheight;
sigma = 2 * sqrt(1/(2*(abs(modelParam(8)))));
MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = sigma * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = residGSH;

% Calculate SNR of GSH signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GSHheight) / noiseSigma_DIFF;

% MM (200728)
% MRS_struct.out.(vox{kk}).(target{jj}).FitError2(ii) = sqrt(mean(residGSH.^2)) / (0.5*noiseSigma_DIFF);

modelFit.modelParam.full = modelParam;
modelFit.freqBounds      = freqBounds;
modelFit.plotBounds      = plotBounds;
modelFit.residPlot       = residPlot;

modelFit.modelParam.Gauss1 = Gauss1;
modelFit.modelParam.Gauss2 = Gauss2;
modelFit.modelParam.Gauss3 = Gauss3;
modelFit.modelParam.Gauss4 = Gauss4;
modelFit.modelParam.Gauss5 = Gauss5;
modelFit.modelParam.Gauss6 = Gauss6;
modelFit.modelParam.Gauss7 = Gauss7;
modelFit.modelParam.Gauss8 = Gauss8;
