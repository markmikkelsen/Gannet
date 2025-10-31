function [MRS_struct, freqBounds, plotBounds, residPlot] = FitGABA(MRS_struct, freq, DIFF, vox, target, ii, jj, kk, lsqopts, nlinopts)

if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target, {'EtOH','GABA','GSH'}))

    freqBounds = find(freq <= 3.55 & freq >= 2.6);
    plotBounds = find(freq <= 3.6 & freq >= 2.6);

    maxinGABA = abs(max(real(DIFF(ii,freqBounds))) - min(real(DIFF(ii,freqBounds))));
    grad_points = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
    LinearInit = grad_points ./ abs(freq(1) - freq(2));
    constInit = (real(DIFF(ii,freqBounds(end))) + real(DIFF(ii,freqBounds(1))))./2;

    GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
    GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow

    lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
    ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
    lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
    ub([1 4 5]) = ub([1 4 5]) / maxinGABA;

    % Down-weight co-edited Cho signal by including
    % observation weights in nonlinear regression
    w = ones(size(DIFF(ii,freqBounds)));
    residfreq = freq(freqBounds);
    ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
    weightRange = ChoRange;
    w(weightRange) = 0.001;

    % Least-squares model fitting
    GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGABA, lb, ub, lsqopts);
    modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
    [GaussModelParam, resid] = nlinfit(freq(freqBounds), sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
    [~, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure

    % Rescale fit parameters and residuals
    GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
    resid = resid * maxinGABA;
    residPlot = residPlot * maxinGABA;

else

    freqBounds = find(freq <= 3.55 & freq >= 2.79);
    plotBounds = find(freq <= 3.6 & freq >= 2.7);

    maxinGABA = abs(max(real(DIFF(ii,freqBounds))) - min(real(DIFF(ii,freqBounds))));
    grad_points = (real(DIFF(ii,freqBounds(end))) - real(DIFF(ii,freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
    LinearInit = grad_points ./ abs(freq(1) - freq(2));
    constInit = (real(DIFF(ii,freqBounds(end))) + real(DIFF(ii,freqBounds(1))))./2;

    GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
    GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow

    lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
    ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
    lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
    ub([1 4 5]) = ub([1 4 5]) / maxinGABA;

    % Down-weight Cho subtraction artifact by including
    % observation weights in nonlinear regression; improves
    % accuracy of peak fitting (MM: 170701 - thanks to Alex
    % Craven of University of Bergen for this idea)
    w = ones(size(DIFF(ii,freqBounds)));
    residfreq = freq(freqBounds);
    ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
    weightRange = ChoRange;
    w(weightRange) = 0.001;

    % Weighted least-squares model fitting
    GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGABA, lb, ub, lsqopts);
    modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
    [GaussModelParam, resid] = nlinfit(freq(freqBounds), sqrt(w) .* real(DIFF(ii,freqBounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
    [~, residPlot] = nlinfit(freq(freqBounds), real(DIFF(ii,freqBounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure

    % Rescale fit parameters and residuals
    GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
    resid = resid * maxinGABA;
    residPlot = residPlot * maxinGABA;

end

GABAheight = GaussModelParam(1);
MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / GABAheight;
MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = GaussModelParam(1) ./ sqrt(-GaussModelParam(2)) * sqrt(pi);
sigma = sqrt(1/(2*(abs(GaussModelParam(2)))));
MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2 * MRS_struct.p.LarmorFreq(ii)) * sigma);
MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;

% Calculate SNR of GABA signal
noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GABAheight) / noiseSigma_DIFF;
