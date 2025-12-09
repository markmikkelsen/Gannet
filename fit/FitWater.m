function [MRS_struct, hb] = FitWater(MRS_struct, freq, WaterData, OFF, vox, target, ii, jj, kk, lsqopts, nlinopts)

if ~isMATLABReleaseOlderThan("R2025a") && MRS_struct.p.append
    font_size_adj  = 2.75;
else
    font_size_adj  = 0;
end

% Estimate height and baseline from data
[maxinWater, watermaxindex] = max(real(WaterData(ii,:)),[],2);
waterbase = mean(real(WaterData(ii,freq <= 4 & freq >= 3.8)));

% Philips data do not phase well based on first point, so do a preliminary
% fit, then adjust phase of WaterData accordingly
% Do this only if the eddy current correction is not performed on water
if strcmp(MRS_struct.p.vendor,'Philips') && ~MRS_struct.p.water_ECC
    % Run preliminary fit of data
    LGModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50];
    freqbounds = freq <= 5.6 & freq >= 3.8;

    % Do the water fit (Lorentz-Gauss)
    LGModelParam = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModel, LGModelInit, nlinopts);

    Eerror = zeros([120 1]);
    for ll = 1:120
        Data = WaterData(ii,freqbounds)*exp(1i*pi/180*ll*3);
        Model = LorentzGaussModel(LGModelParam,freq(freqbounds));
        Eerror(ll) = sum((real(Data)-Model).^2);
    end
    [~,index] = min(Eerror);
    WaterData(ii,:) = WaterData(ii,:) * exp(1i*pi/180*index*3);
end

LGPModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50 0];
lb = [0.01*maxinWater 1   4.6 0     0 -50 -pi];
ub = [40*maxinWater   100 4.8 1e-6  1  0   pi];

freqbounds = freq <= 5.6 & freq >= 3.8;

% Least-squares model fitting
LGPModelInit = lsqcurvefit(@LorentzGaussModelP, LGPModelInit, freq(freqbounds), real(WaterData(ii,freqbounds)), lb, ub, lsqopts);
[LGPModelParam, residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModelP, LGPModelInit, nlinopts);

WaterArea = sum(real(LorentzGaussModel(LGPModelParam(1:end-1), freq(freqbounds))) - BaselineModel(LGPModelParam(3:5), freq(freqbounds)),2);
MRS_struct.out.(vox{kk}).water.Area(ii) = WaterArea * abs(freq(1) - freq(2));
waterheight = LGPModelParam(1);
MRS_struct.out.(vox{kk}).water.FitError(ii) = 100 * std(residw) / waterheight;

LG = real(LorentzGaussModel(LGPModelParam(1:end-1), freq(freqbounds))) - BaselineModel(LGPModelParam(3:5), freq(freqbounds));
LG = LG/max(LG);
ind = find(LG >= 0.5);
f = freq(freqbounds);
w = abs(f(ind(1)) - f(ind(end)));
MRS_struct.out.(vox{kk}).water.FWHM(ii) = w * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).water.ModelParam(ii,:) = LGPModelParam;
MRS_struct.out.(vox{kk}).water.Resid(ii,:) = residw;

% Calculate SNR of water signal
noiseSigma_Water = CalcNoise(freq, WaterData(ii,:));
MRS_struct.out.(vox{kk}).water.SNR(ii) = abs(waterheight) / noiseSigma_Water;

% Water spectrum plot
hb = subplot(2,2,3);
watmin = min(real(WaterData(ii,:)));
watmax = max(real(WaterData(ii,:)));
residw = residw + watmin - max(residw);
plot(freq(freqbounds), real(WaterData(ii,freqbounds)), 'b', ...
    freq(freqbounds), real(LorentzGaussModelP(LGPModelParam,freq(freqbounds))), 'r', ...
    freq(freqbounds), residw, 'k');
set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XTick', 4.2:0.2:5.2, 'FontSize', 10 - font_size_adj);
xlim([4.2 5.2]);
set(get(gca,'YAxis'),'Visible','off');
% Add on some labels
text(4.8, watmax/2, 'Water', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
labelfreq = freq(freqbounds);
rlabelbounds = labelfreq <= 4.4 & labelfreq >= 4.25;
axis_bottom = axis;
text(4.4, max(min(residw(rlabelbounds)) - 0.05 * watmax, axis_bottom(3)), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
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

    residWater = OFF(ii,:);
    freqWater  = freq <= 4.68 + 0.4 & freq >= 4.68 - 0.4;
    residWater = residWater(freqWater);
    freqWater  = freq(freqWater);

    [maxResidWater, maxInd] = max(abs(real(residWater)));
    s = sign(real(residWater(maxInd)));
    maxResidWater = s * maxResidWater;
    offset = real(residWater(1));

    LGPModelInit = [maxResidWater 25 freqWater(maxInd) 0 offset 0.001 0];
    LGPModelInit([1 5]) = LGPModelInit([1 5]) / maxResidWater; % Scale initial conditions to avoid warnings about numerical underflow

    %lb = [maxResidWater-abs(2*maxResidWater) 1 freqWaterOFF(maxInd)-0.2 0 0 -200 -pi];
    %ub = [maxResidWater+abs(2*maxResidWater) 100 freqWaterOFF(maxInd)+0.2 0.000001 1 0 pi];
    %lb([1 4 5]) = lb([1 4 5]) / maxResidWater;
    %ub([1 4 5]) = ub([1 4 5]) / maxResidWater;

    % Least-squares model fitting
    [MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,:), residRW] = nlinfit(freqWater, real(residWater) / maxResidWater, @LorentzGaussModelP, LGPModelInit, nlinopts);

    % Rescale fit parameters and residuals
    MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) = MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) * maxResidWater;
    residRW = residRW * maxResidWater;

    MRS_struct.out.(vox{kk}).resid_water.FitError(ii) = 100*std(residRW) / MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1);

    MRS_struct.out.(vox{kk}).resid_water.SuppressionFactor(ii) = ...
        100 * ((MRS_struct.out.(vox{kk}).water.ModelParam(ii,1) - abs(MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1))) ...
        / MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));

end
