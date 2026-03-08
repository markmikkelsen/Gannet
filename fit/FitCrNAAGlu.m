function [MRS_struct, freqBoundsChoCr, freqBoundsCr, residCr] = ...
    FitCrNAAGlu(MRS_struct, freq, OFF, SUM, vox, ii, kk, baseline, lsqopts, nlinopts, lsqnlinopts)

%% Cr Fit

freqBoundsChoCr = find(freq <= 3.6 & freq >= 2.6);
Baseline_offset = real(OFF(freqBoundsChoCr(1)) + OFF(freqBoundsChoCr(end))) / 2;
Width_estimate  = 0.03;
Area_estimate   = (max(real(OFF(freqBoundsChoCr))) - min(real(OFF(freqBoundsChoCr)))) * Width_estimate * 4;
ChoCrModelParamInit = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1] ...
                        .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
[ChoCrModelParam, ~, residChoCr] = FitChoCr(freq(freqBoundsChoCr), OFF(freqBoundsChoCr), ChoCrModelParamInit, MRS_struct.p.LarmorFreq(ii));
MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:) = ChoCrModelParam ./ [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];

% Initialise fitting pars
freqBoundsCr = find(freq <= 3.12 & freq >= 2.72);
CrModelParamInit = [max(real(OFF(freqBoundsCr))) 0.05 3.0 0 0 0];

% Least-squares model fitting
CrModelParamInit = lsqcurvefit(@LorentzModel, CrModelParamInit, freq(freqBoundsCr), real(OFF(freqBoundsCr)), [], [], lsqopts);
[CrModelParam, residCr] = nlinfit(freq(freqBoundsCr), real(OFF(freqBoundsCr)), @LorentzModel, CrModelParamInit, nlinopts);

MRS_struct.out.(vox{kk}).Cr.ModelParam(ii,:) = CrModelParam;
CrHeight = CrModelParam(1) / (2*pi*CrModelParam(2));
MRS_struct.out.(vox{kk}).Cr.FitError(ii)     = 100 * std(residCr) / CrHeight;
MRS_struct.out.(vox{kk}).Cr.Resid(ii,:)      = residCr;

MRS_struct.out.(vox{kk}).Cho.ModelParam(ii,:) = MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:);
ChoHeight = (MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1) * MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,7)) / ...
                (2 * pi * MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2));
MRS_struct.out.(vox{kk}).Cho.FitError(ii) = 100 * std(residChoCr) / ChoHeight;
MRS_struct.out.(vox{kk}).Cho.Resid(ii,:) = residChoCr;

MRS_struct.out.(vox{kk}).Cr.Area(ii)  = sum(real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqBoundsChoCr)) - ...
    TwoLorentzModel([0 MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2:end-1) 0], freq(freqBoundsChoCr)))) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).Cho.Area(ii) = sum(real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:), freq(freqBoundsChoCr)) - ...
    TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqBoundsChoCr)))) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).Cr.FWHM(ii)  = ChoCrModelParam(2);
MRS_struct.out.(vox{kk}).Cho.FWHM(ii) = ChoCrModelParam(2);

% Calculate SNR of Cr signal
noiseSigma_OFF = CalcNoise(freq, OFF);
MRS_struct.out.(vox{kk}).Cr.SNR(ii)  = abs(CrHeight) / noiseSigma_OFF;
MRS_struct.out.(vox{kk}).Cho.SNR(ii) = abs(ChoHeight) / noiseSigma_OFF;


%% NAA Fit

freqBoundsNAA = find(freq <= 2.25 & freq >= 1.75);
maxinNAA   = max(real(OFF(freqBoundsNAA)));
gradPoints = (real(OFF(freqBoundsNAA(end))) - real(OFF(freqBoundsNAA(1)))) ./ abs(freqBoundsNAA(end) - freqBoundsNAA(1));
LinearInit = gradPoints ./ abs(freq(1) - freq(2));
constInit  = (real(OFF(freqBoundsNAA(end))) + real(OFF(freqBoundsNAA(1)))) ./ 2;

NAAModelParamInit = [maxinNAA 0.05 2.01 0 -LinearInit constInit];
lb = [0 0.01 1.97 0 -40*maxinNAA -2000*maxinNAA];
ub = [4000*maxinNAA 0.1 2.05 0.5 40*maxinNAA 1000*maxinNAA];

% Least-squares model fitting
NAAModelParamInit = lsqcurvefit(@LorentzModel, NAAModelParamInit, freq(freqBoundsNAA), real(OFF(freqBoundsNAA)), lb, ub, lsqopts);
[NAAModelParam, resid] = nlinfit(freq(freqBoundsNAA), real(OFF(freqBoundsNAA)), @LorentzModel, NAAModelParamInit, nlinopts);

NAAheight = NAAModelParam(1) / (2 * pi * NAAModelParam(2));
MRS_struct.out.(vox{kk}).NAA.FitError(ii) = 100 * std(resid) / NAAheight;
NAAModelParam_tmp = NAAModelParam;
NAAModelParam_tmp(4:6) = 0;
MRS_struct.out.(vox{kk}).NAA.Area(ii) = sum(LorentzModel(NAAModelParam_tmp, freq(freqBoundsNAA))) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).NAA.FWHM(ii) = abs((2 * MRS_struct.p.LarmorFreq(ii)) * NAAModelParam_tmp(2));
MRS_struct.out.(vox{kk}).NAA.ModelParam(ii,:) = NAAModelParam;
MRS_struct.out.(vox{kk}).NAA.Resid(ii,:) = resid;

% Calculate SNR of NAA signal
MRS_struct.out.(vox{kk}).NAA.SNR(ii) = abs(NAAheight) / noiseSigma_OFF;


%% Glu Fit
% 2025-11: now fit Glu in the SUM; this still needs testing to see if it's reliable
% 2025-12: new fitting approach using predefined baseline

% Glu_fitLim = [2.34-0.1 2.34+0.1];
Glu_fitLim = [2.22 2.47];
freqBounds = find(freq <= Glu_fitLim(2) & freq >= Glu_fitLim(1));

maxinGlu = max(real(SUM(freqBounds)));
widthGlu = 30;

GluModelParamInit = [maxinGlu 0.5 0.5 ... % amplitudes
                    widthGlu ... % width
                    2.34 0.04 ... % freqs / J-coupling constant
                    0]; % phase

lb = [0 0 0 ...
      0 ...
      2.34-0.02 0 ...
      0];
ub = [2*maxinGlu 2 2 ...
      150 ...
      2.34+0.02 0.05 ...
      pi];

% Scale initial conditions and bounds to avoid warnings about numerical underflow
amplParams = [1 2 3];
GluModelParamInit(amplParams(1)) = GluModelParamInit(amplParams(1)) / maxinGlu;
lb(amplParams(1)) = lb(amplParams(1)) / maxinGlu;
ub(amplParams(1)) = ub(amplParams(1)) / maxinGlu;

% Model fitting with fixed baseline
[GluModelParam, resid, h_tmp] = FitSignalModel(@ThreeLorentzModel_phased_noBaseline, ... % model
                                freq(freqBounds), ... % freq
                                real(SUM(freqBounds)) / maxinGlu, ... % data
                                baseline(freqBounds) / maxinGlu, ... % baseline
                                GluModelParamInit, ... % beta0
                                lb, ...
                                ub, ...
                                lsqnlinopts);

% Rescale fit parameters and residuals
GluModelParam(amplParams(1)) = GluModelParam(amplParams(1)) * maxinGlu;
resid = resid * maxinGlu;

% Lorentzians
[Lorentz1, Lorentz2, Lorentz3] = deal(GluModelParam);
Lorentz1(amplParams(2:3)) = 0;
Lorentz2(amplParams(2))   = 0;
Lorentz3(amplParams(3))   = 0;

% Plot individual Lorentzians
hold on;
plot(freq, real(SUM) / maxinGlu, 'k');
plot(freq(freqBounds), ThreeLorentzModel_phased_noBaseline(Lorentz1, freq(freqBounds)) / maxinGlu);
plot(freq(freqBounds), (ThreeLorentzModel_phased_noBaseline(Lorentz2, freq(freqBounds)) - ...
    ThreeLorentzModel_phased_noBaseline(Lorentz1, freq(freqBounds))) / maxinGlu);
plot(freq(freqBounds), (ThreeLorentzModel_phased_noBaseline(Lorentz3, freq(freqBounds)) - ...
    ThreeLorentzModel_phased_noBaseline(Lorentz1, freq(freqBounds))) / maxinGlu);
legend({'data','model + baseline','baseline','residual', ...
    'lorentz1','lorentz2','lorentz3'}, ...
    'Box','off','Location','best');
xlim([1.75 3.5]);
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
% exportgraphics(h_tmp, fullfile(out_dir, [fname '_Glu_model_fit.png']), "Resolution", 300);
close(h_tmp);

GluHeight = GluModelParam(1) * GluModelParam(4);
MRS_struct.out.(vox{kk}).Glu.FitError(ii)     = 100 * std(resid) / GluHeight;
MRS_struct.out.(vox{kk}).Glu.Area(ii)         = GluModelParam(1) * pi;
MRS_struct.out.(vox{kk}).Glu.FWHM(ii)         = 2 / GluModelParam(4) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).Glu.ModelParam(ii,:) = GluModelParam;
MRS_struct.out.(vox{kk}).Glu.Resid(ii,:)      = resid;

% Calculate SNR of Glu signal
noiseSigma_SUM = CalcNoise(freq, SUM);
MRS_struct.out.(vox{kk}).Glu.SNR(ii) = abs(GluHeight) / noiseSigma_SUM;
