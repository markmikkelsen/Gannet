function [MRS_struct, hb] = FitWater(MRS_struct, freq, water, OFF, vox, target, ii, jj, kk, lsqopts, nlinopts)

if ~isMATLABReleaseOlderThan("R2025a") && MRS_struct.p.append
    font_size_adj  = 2.75;
else
    font_size_adj  = 0;
end

% Estimate height and baseline from data
[maxinWater, waterMaxInd] = max(real(water));
offset = mean(real(water(freq <= 4 & freq >= 3.8)));

% Philips data do not phase well based on first point, so do a preliminary
% fit, then adjust phase of water accordingly
% Do this only if the eddy current correction is not performed on water
if strcmp(MRS_struct.p.vendor,'Philips') && ~MRS_struct.p.water_ECC
    % Run preliminary fit of data
    modelParamInit = [maxinWater 20 freq(waterMaxInd) 0 offset -50];
    freqBounds = freq <= 5.6 & freq >= 3.8;

    % Do the water fit (Lorentz-Gauss)
    modelParam = nlinfit(freq(freqBounds), real(water(freqBounds)), @LorentzGaussModel, modelParamInit, nlinopts);

    Eerror = zeros([120 1]);
    for ll = 1:120
        Data = water(freqBounds) * exp(1i * pi/180 * ll * 3);
        Model = LorentzGaussModel(modelParam,freq(freqBounds));
        Eerror(ll) = sum((real(Data) - Model).^2);
    end
    [~,index] = min(Eerror);
    water = water * exp(1i * pi/180 * index * 3);
end

modelParamInit = [maxinWater 20 freq(waterMaxInd) 0 offset -20 0];
lb = [0.01*maxinWater 1   4.6 0    0 -100 -pi];
ub = [40*maxinWater   100 4.8 1e-6 1  0    pi];

freqBounds = freq <= 5.6 & freq >= 3.8;

% Least-squares model fitting
modelParamInit = lsqcurvefit(@LorentzGaussModelP, modelParamInit, freq(freqBounds), real(water(freqBounds)), lb, ub, lsqopts);
[modelParam, resid] = nlinfit(freq(freqBounds), real(water(freqBounds)), @LorentzGaussModelP, modelParamInit, nlinopts);

waterModelParam = modelParam(1:end-1); % exclude phase parameter
waterModelParam(4:5) = 0; % zero baseline

MRS_struct.out.(vox{kk}).water.Area(ii) = sum(LorentzGaussModel(waterModelParam, freq(freqBounds))) * abs(freq(1) - freq(2));
waterHeight = modelParam(1);
MRS_struct.out.(vox{kk}).water.FitError(ii) = 100 * std(resid) / waterHeight;

LG = real(LorentzGaussModel(waterModelParam, freq(freqBounds)));
LG = LG / max(LG);
ind = find(LG >= 0.5);
f = freq(freqBounds);
w = abs(f(ind(1)) - f(ind(end)));
MRS_struct.out.(vox{kk}).water.FWHM(ii) = w * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).water.ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).water.Resid(ii,:) = resid;

% Calculate SNR of water signal
noiseSigma_Water = CalcNoise(freq, water);
MRS_struct.out.(vox{kk}).water.SNR(ii) = abs(waterHeight) / noiseSigma_Water;

% Water spectrum plot
hb = subplot(2,2,3);
watmin = min(real(water));
watmax = max(real(water));
resid = resid + watmin - max(resid);
plot(freq(freqBounds), real(water(freqBounds)), 'b', ...
    freq(freqBounds), real(LorentzGaussModelP(modelParam, freq(freqBounds))), 'r', ...
    freq(freqBounds), resid, 'k');
set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XTick', 4.2:0.2:5.2, 'FontSize', 10 - font_size_adj);
xlim([4.2 5.2]);
set(get(gca,'YAxis'),'Visible','off');
% Add on some labels
text(4.8, watmax/2, 'Water', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
labelfreq = freq(freqBounds);
rlabelbounds = labelfreq <= 4.4 & labelfreq >= 4.25;
axis_bottom = axis;
text(4.4, max(min(resid(rlabelbounds)) - 0.05 * watmax, axis_bottom(3)), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
xlabel('ppm', 'FontSize', 11 - font_size_adj);
title('Reference signals', 'FontSize', 11 - font_size_adj);

% Root sum square fit error and concentration in institutional units
switch target{jj}
    case 'GABA'
        MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);

    case 'Glx'
        MRS_struct.out.(vox{kk}).Glx.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);

    case 'GABAGlx'
        MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct.out.(vox{kk}).Glx.FitError_W(ii)  = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);

    case 'GSH'
        MRS_struct.out.(vox{kk}).GSH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GSH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);

    case 'Lac'
        MRS_struct.out.(vox{kk}).Lac.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Lac.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);

    case 'EtOH'
        MRS_struct.out.(vox{kk}).EtOH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).EtOH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
        MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
end

% Generate scaled spectra (for plotting)
MRS_struct.spec.(vox{kk}).(target{jj}).off_scaled(ii,:) = ...
    MRS_struct.spec.(vox{kk}).(target{jj}).off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
if MRS_struct.p.HERMES
    MRS_struct.spec.(vox{kk}).(target{jj}).off_off_scaled(ii,:) = ...
        MRS_struct.spec.(vox{kk}).(target{jj}).off_off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
end
MRS_struct.spec.(vox{kk}).(target{jj}).on_scaled(ii,:) = ...
    MRS_struct.spec.(vox{kk}).(target{jj}).on(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
MRS_struct.spec.(vox{kk}).(target{jj}).diff_unfilt_scaled(ii,:) = ...
    MRS_struct.spec.(vox{kk}).(target{jj}).diff_unfilt(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
MRS_struct.spec.(vox{kk}).(target{jj}).diff_scaled(ii,:) = ...
    MRS_struct.spec.(vox{kk}).(target{jj}).diff(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
MRS_struct.spec.(vox{kk}).(target{jj}).sum_scaled(ii,:) = ...
    MRS_struct.spec.(vox{kk}).(target{jj}).sum(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));

% Reorder structure fields
MRS_struct.out.(vox{kk}).water = orderfields(MRS_struct.out.(vox{kk}).water, {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2a.  Residual Water Fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if MRS_struct.p.fit_resid_water

    freqBounds = freq <= 4.68 + 0.4 & freq >= 4.68 - 0.4;
    residWater = OFF(freqBounds);
    freqWater  = freq(freqBounds);

    [maxResidWater, maxInd] = max(abs(real(residWater)));
    s = sign(real(residWater(maxInd)));
    maxResidWater = s * maxResidWater;
    offset = real(residWater(1));

    modelParamInit = [maxResidWater 25 freqWater(maxInd) 0 offset 0.001 0];
    modelParamInit([1 5]) = modelParamInit([1 5]) / maxResidWater; % scale initial conditions to avoid warnings about numerical underflow

    % Nonlinear least-squares model fitting
    [MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,:), resid] = ...
        nlinfit(freqWater, real(residWater) / maxResidWater, @LorentzGaussModelP, modelParamInit, nlinopts);

    % Rescale fit parameters and residuals
    MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) = MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) * maxResidWater;
    resid = resid * maxResidWater;

    MRS_struct.out.(vox{kk}).resid_water.FitError(ii) = 100 * std(resid) / MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1);

    MRS_struct.out.(vox{kk}).resid_water.SuppressionFactor(ii) = ...
        100 * ((MRS_struct.out.(vox{kk}).water.ModelParam(ii,1) - abs(MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1))) / ...
        MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));

end
