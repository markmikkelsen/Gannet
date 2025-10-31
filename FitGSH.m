function [MRS_struct, modelFit] = FitGSH(MRS_struct, freq, DIFF, vox, target, ii, jj, kk, lsqopts, nlinopts)

freqBounds = find(freq <= 3.5 & freq >= 2.25);
plotBounds = find(freq <= 4.2 & freq >= 1.75);

GSHbounds = freq <= 3.3 & freq >= 2.85;
AspBounds = freq <= 2.85 & freq >= 2.25;

maxinGSH = max(abs(real(DIFF(ii,GSHbounds))));
[maxinAsp, maxInd] = max(abs(real(DIFF(ii,AspBounds))));

offset     = real(DIFF(ii,freqBounds(end)));
gradPoints = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
linearInit = gradPoints ./ abs(freq(1) - freq(2));

DIFF_Asp = DIFF(ii,AspBounds);
s = sign(real(DIFF_Asp(maxInd)));
maxinAsp = s * maxinAsp;

if MRS_struct.p.HERMES
    s = -1;
else
    s = 1;
end

if MRS_struct.p.TE(ii) < 100

        GSHgaussModel = @EightGaussModel;

        GaussModelInit = [maxinGSH        -300  2.95 ...
                         -maxinGSH*0.1    -300  2.77 ...            
                          s*maxinAsp*0.25 -500  2.73 ...
                          maxinAsp        -1000 2.63 ...
                          maxinAsp        -1000 2.57 ...
                          maxinAsp        -1000 2.50 ...
                          s*maxinAsp*0.15 -600  2.45 ...
                         -maxinGSH*0.1 -300  2.35 ...
                          offset -linearInit -linearInit];

        lb = [-4000*maxinGSH       -1000 2.95-0.02 ...
              -4000*maxinGSH*0.1   -1000 2.77-0.02 ...            
               4000*maxinAsp*0.25  -1000 2.73-0.02 ...
               4000*maxinAsp       -1000 2.63-0.02 ...
               4000*maxinAsp       -1000 2.57-0.02 ...
               4000*maxinAsp       -1000 2.50-0.02 ...
               4000*maxinAsp*0.15  -1000 2.45-0.02 ...
              -4000*maxinGSH*0.1   -1000 2.35-0.02 ...
              -2000*abs(offset) 2000*maxinAsp 2000*maxinAsp];
        ub =  [4000*maxinGSH      -40 2.95+0.02 ...
               4000*maxinGSH*0.1  -40 2.77+0.02 ...            
              -4000*maxinAsp*0.25 -40 2.73+0.02 ...
              -4000*maxinAsp      -40 2.63+0.02 ...
              -4000*maxinAsp      -40 2.57+0.02 ...
              -4000*maxinAsp      -40 2.50+0.02 ...
              -4000*maxinAsp*0.15 -40 2.45+0.02 ...
               4000*maxinGSH*0.1  -40 2.35+0.02 ...        
               1000*abs(offset) -1000*maxinAsp -1000*maxinAsp];

else

    GSHgaussModel = @SevenGaussModel;

    GaussModelInit = [maxinGSH        -300  2.95 ...
                     -maxinGSH*0.25   -1000 2.78 ...
                      s*maxinAsp*0.25 -300  2.72 ...
                      maxinAsp        -1000 2.63 ...
                      maxinAsp        -1000 2.57 ...
                      maxinAsp        -1000 2.45 ...
                      s*maxinAsp*0.15 -1000 2.38 ...
                      offset -linearInit -linearInit];

    lb = [-4*maxinGSH      -1000 2.95-0.02 ...
          -4*maxinGSH*0.25 -4000 2.78-0.02 ...
           4*maxinAsp*0.25 -1000 2.72-0.02 ...
           4*maxinAsp      -1000 2.63-0.02 ...
           4*maxinAsp      -1000 2.57-0.02 ...
           4*maxinAsp      -1000 2.45-0.02 ...
           4*maxinAsp*0.15 -4000 2.38-0.02 ...
          -10*abs(offset) 10*maxinAsp 10*maxinAsp];
    ub =  [4*maxinGSH      -40 2.95+0.02 ...
           4*maxinGSH*0.25 -40 2.78+0.02 ...
          -4*maxinAsp*0.25 -40 2.72+0.02 ...
          -4*maxinAsp      -40 2.63+0.02 ...
          -4*maxinAsp      -40 2.57+0.02 ...
          -4*maxinAsp      -40 2.45+0.02 ...
          -4*maxinAsp*0.15 -40 2.38+0.02 ...
           10*abs(offset) -10*maxinAsp -10*maxinAsp];

end

% Scale initial conditions and bounds to avoid warnings about numerical underflow
if MRS_struct.p.TE(ii) < 100
    GaussModelInit([1 4 7 10 13 16 19 22 25 26 27]) = ... % eight-Gaussian model
        GaussModelInit([1 4 7 10 13 16 19 22 25 26 27]) / maxinGSH;
    lb([1 4 7 10 13 16 19 22 25 26 27]) = ...
        lb([1 4 7 10 13 16 19 22 25 26 27]) / maxinGSH;
    ub([1 4 7 10 13 16 19 22 25 26 27]) = ...
        ub([1 4 7 10 13 16 19 22 25 26 27]) / maxinGSH;
else
    GaussModelInit([1 4 7 10 13 16 19 22 23 24]) = ... % seven-Gaussian model
        GaussModelInit([1 4 7 10 13 16 19 22 23 24]) / maxinGSH;
    lb([1 4 7 10 13 16 19 22 23 24]) = ...
        lb([1 4 7 10 13 16 19 22 23 24]) / maxinGSH;
    ub([1 4 7 10 13 16 19 22 23 24]) = ...
        ub([1 4 7 10 13 16 19 22 23 24]) / maxinGSH;
end

if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target, {'EtOH','GABA','GSH'}))
    w = ones(size(DIFF(ii,freqBounds)));
    residfreq = freq(freqBounds);
    ChoRange = residfreq >= 3.13 & residfreq <= 3.3;
    weightRange = ChoRange;
    w(weightRange) = 0.001;
else
    w = ones(size(DIFF(ii,freqBounds)));
end

% Weighted least-squares model fitting
GaussModelInit = lsqcurvefit(GSHgaussModel, GaussModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGSH, lb, ub, lsqopts);
modelFun_w = @(x,freq) sqrt(w) .* GSHgaussModel(x,freq); % add weights to the model
[modelParam, resid] = nlinfit(freq(freqBounds), sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinGSH, ...
                                modelFun_w, GaussModelInit, nlinopts); % add weights to the data
[~, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGSH, ...
                    GSHgaussModel, modelParam, nlinopts); % re-run for residuals for output figure

% Rescale fit parameters and residuals
if MRS_struct.p.TE(ii) < 100
    modelParam([1 4 7 10 13 16 19 22 25 26 27]) = modelParam([1 4 7 10 13 16 19 22 25 26 27]) * maxinGSH;
else
    modelParam([1 4 7 10 13 16 19 22 23 24]) = modelParam([1 4 7 10 13 16 19 22 23 24]) * maxinGSH;
end
resid = resid * maxinGSH;
residPlot = residPlot * maxinGSH;

% %%%
% GSHgaussModelParam = modelParam;
% GSHgaussModelParam([4:3:22 25:27]) = 0;
% figure;
% hold on;
% PlotSpec(freq, real(DIFF(ii,:)));
% PlotSpec(freq(freqBounds), GSHgaussModel(modelParam, freq(freqBounds)));
% hold off;
% yyaxis right;
% PlotSpec(freq(freqBounds), GSHgaussModel(GSHgaussModelParam, freq(freqBounds)));
% set(gca,'XLim',[2.25 3.25]);
% %%%

% Range to determine residuals for GSH
residfreq = freq(freqBounds);
residGSH  = resid(residfreq <= 3.3 & residfreq >= 2.82);

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

modelFit.modelParam  = modelParam;
modelFit.freqBounds  = freqBounds;
modelFit.plotBounds  = plotBounds;
modelFit.residPlot   = residPlot;
if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target, {'EtOH','GABA','GSH'}))
    modelFit.weightRange = weightRange;
end
