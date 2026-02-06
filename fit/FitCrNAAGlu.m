function [MRS_struct, Cr_OFF, freqBoundsChoCr, freqBoundsCr, residCr] = ...
    FitCrNAAGlu(MRS_struct, freq, OFF, SUM, vox, ii, kk, baseline, lsqopts, nlinopts, lsqnlinopts)
                
% Cr Fit

Cr_OFF = OFF(ii,:);
freqBoundsChoCr = freq <= 3.6 & freq >= 2.6;

Cho_Cr = Cr_OFF(freqBoundsChoCr).';
Baseline_offset = real(Cho_Cr(1) + Cho_Cr(end)) / 2;
Width_estimate  = 0.05;
Area_estimate   = (max(real(Cho_Cr)) - min(real(Cho_Cr))) * Width_estimate * 4;
ChoCrModelInit  = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1] ...
                    .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
[ChoCrModelParam, ~, residChoCr] = FitChoCr(freq(freqBoundsChoCr), Cho_Cr, ChoCrModelInit, MRS_struct.p.LarmorFreq(ii));
MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:) = ChoCrModelParam ./ [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];

% Initialise fitting pars
freqBoundsCr = freq <= 3.12 & freq >= 2.72;
modelParamInit = [max(real(Cr_OFF(freqBoundsCr))) 0.05 3.0 0 0 0];

% Least-squares model fitting
modelParamInit = lsqcurvefit(@LorentzModel, modelParamInit, freq(freqBoundsCr), real(Cr_OFF(freqBoundsCr)), [], [], lsqopts);
[modelParam, residCr] = nlinfit(freq(freqBoundsCr), real(Cr_OFF(freqBoundsCr)), @LorentzModel, modelParamInit, nlinopts);

MRS_struct.out.(vox{kk}).Cr.ModelParam(ii,:) = modelParam;
CrHeight = modelParam(1) / (2*pi*modelParam(2));
MRS_struct.out.(vox{kk}).Cr.FitError(ii)     = 100 * std(residCr) / CrHeight;
MRS_struct.out.(vox{kk}).Cr.Resid(ii,:)      = residCr;

MRS_struct.out.(vox{kk}).Cho.ModelParam(ii,:) = MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:);
ChoHeight = (MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1) * MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,7)) / (2*pi*MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2));
MRS_struct.out.(vox{kk}).Cho.FitError(ii)     = 100 * std(residChoCr) / ChoHeight;
MRS_struct.out.(vox{kk}).Cho.Resid(ii,:)      = residChoCr;

MRS_struct.out.(vox{kk}).Cr.Area(ii)  = sum(real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqBoundsChoCr)) - ...
    TwoLorentzModel([0 MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2:end-1) 0], freq(freqBoundsChoCr)))) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).Cho.Area(ii) = sum(real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:), freq(freqBoundsChoCr)) - ...
    TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqBoundsChoCr)))) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).Cr.FWHM(ii)  = ChoCrModelParam(2);
MRS_struct.out.(vox{kk}).Cho.FWHM(ii) = ChoCrModelParam(2);

% Calculate SNR of Cr signal
noiseSigma_OFF = CalcNoise(freq, OFF(ii,:));
MRS_struct.out.(vox{kk}).Cr.SNR(ii)  = abs(CrHeight) / noiseSigma_OFF;
MRS_struct.out.(vox{kk}).Cho.SNR(ii) = abs(ChoHeight) / noiseSigma_OFF;


% NAA Fit

NAA_OFF = OFF(ii,:);
freqBounds = find(freq <= 2.25 & freq >= 1.75);

maxinNAA   = max(real(NAA_OFF(freqBounds)));
gradPoints = (real(NAA_OFF(freqBounds(end))) - real(NAA_OFF(freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
LinearInit = gradPoints ./ abs(freq(1) - freq(2));
constInit  = (real(NAA_OFF(freqBounds(end))) + real(NAA_OFF(freqBounds(1)))) ./ 2;

modelParamInit = [maxinNAA 0.05 2.01 0 -LinearInit constInit];
lb = [0 0.01 1.97 0 -40*maxinNAA -2000*maxinNAA];
ub = [4000*maxinNAA 0.1 2.05 0.5 40*maxinNAA 1000*maxinNAA];

% Least-squares model fitting
modelParamInit = lsqcurvefit(@LorentzModel, modelParamInit, freq(freqBounds), real(NAA_OFF(freqBounds)), lb, ub, lsqopts);
[modelParam, resid] = nlinfit(freq(freqBounds), real(NAA_OFF(freqBounds)), @LorentzModel, modelParamInit, nlinopts);

NAAheight = modelParam(1) / (2*pi*modelParam(2));
MRS_struct.out.(vox{kk}).NAA.FitError(ii) = 100 * std(resid) / NAAheight;
NAAModelParam = modelParam;
NAAModelParam(4) = 0;
MRS_struct.out.(vox{kk}).NAA.Area(ii) = sum(LorentzModel(NAAModelParam,freq(freqBounds)) - BaselineModel(NAAModelParam([3 6 5]),freq(freqBounds)), 2) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).NAA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * NAAModelParam(2));
MRS_struct.out.(vox{kk}).NAA.ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).NAA.Resid(ii,:) = resid;

% Calculate SNR of NAA signal
MRS_struct.out.(vox{kk}).NAA.SNR(ii) = abs(NAAheight) / noiseSigma_OFF;


% Glu Fit
% 2025-11: now fit Glu in the SUM; this still needs testing to see if it's reliable
% 2025-12: new fitting approach using predefined baseline

Glu_SUM = SUM(ii,:);
Glu_fitLim = [2.34-0.04 2.34+0.04];
freqBounds = find(freq <= Glu_fitLim(2) & freq >= Glu_fitLim(1));

maxinGlu = max(real(Glu_SUM(freqBounds)));
% gradPoints = (real(Glu_SUM(freqBounds(end))) - real(Glu_SUM(freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1));
% LinearInit = gradPoints ./ abs(freq(1) - freq(2));
% constInit  = (real(Glu_SUM(freqBounds(end))) + real(Glu_SUM(freqBounds(1)))) ./ 2;
widthGlu = 30;

% modelParamInit = [maxinGlu -0.2 -0.2 widthGlu 2.34 0.05 0];
modelParamInit = [maxinGlu widthGlu 2.34 0];

% lb = [0 -1 -1 0 2.34-0.02 0 0];
% ub = [2*maxinGlu 0 0 200 2.34+0.02 0.1 pi];
lb = [0 0 2.34-0.02 0];
ub = [2*maxinGlu 200 2.34+0.02 pi];

% GaussModel = @(x,freq) x(1) * exp(x(2) * (freq - x(3)).^2);

% GaussModelInit = [ maxinGlu -maxinGlu*0.15 -maxinGlu*0.15 ...
%                   -300 -300 -300 ...
%                    2.34 2.34-0.05 2.34+0.05];
% 
% lb = [ 0 -2*maxinGlu*0.15 -2*maxinGlu*0.15 ...
%       -5000 -5000 -5000 ...
%        2.34-0.02 2.34-0.05-0.02 2.34+0.05-0.02];
% 
% ub = [2*maxinGlu 2*maxinGlu 2*maxinGlu ...
%       -100 -100 -100 ...
%       2.34+0.02 2.34-0.05+0.02 2.34+0.05+0.02];

% GaussModelInit = [maxinGlu -300 2.34];

% lb = [0 -5000 2.34-0.02];
% ub = [2*maxinGlu -500 2.34+0.02];

% a = maxinGlu/widthInit;
% amplParams = [1 4 7];
% LorentzModelInit(amplParams) = LorentzModelInit(amplParams) / a;
% lb(amplParams) = lb(amplParams) / a;
% ub(amplParams) = ub(amplParams) / a;

% Least-squares model fitting
% LorentzModelInit = lsqcurvefit(@ThreeLorentzModel_linBaseline, LorentzModelInit, freq(freqBounds), real(Glu_SUM(freqBounds)), lb, ub, lsqopts);
% [LorentzModelParam, resid] = nlinfit(freq(freqBounds), real(Glu_SUM(freqBounds)), @ThreeLorentzModel_linBaseline, LorentzModelInit, nlinopts);

% Model fitting with predetermined baseline
% baseFreq = freq(freq < 6 & freq > 0);
% baseLim = baseFreq <= Glu_fitLim(2) & baseFreq >= Glu_fitLim(1);
[modelParam, resid, h_tmp] = FitSignalModel(@LorentzModel_phased_noBaseline, ... % model
                                freq(freqBounds), ... % freq
                                real(Glu_SUM(freqBounds)) / maxinGlu, ... % data
                                baseline(freqBounds) / maxinGlu, ... % baseline
                                modelParamInit, ... % beta0
                                lb, ...
                                ub, ...
                                lsqnlinopts);
% [GaussModelParam, resid, h_tmp] = FitSignalModel(GaussModel, ... % model
%                                     freq(freqBounds), ... % freq
%                                     real(Glu_SUM(freqBounds)), ... % data
%                                     baseline(freqBounds), ... % baseline
%                                     GaussModelInit, ... % beta0
%                                     lb, ...
%                                     ub, ...
%                                     lsqnlinopts);

% Rescale fit parameters and residuals
modelParam(1) = modelParam(1) * maxinGlu;
resid = resid * maxinGlu;

% Plot model without baseline
L = LorentzModel_phased_noBaseline(modelParam, freq(freqBounds)) / maxinGlu;

% amplParams = [1 2 3];
% [Gauss1, Gauss2, Gauss3] = deal(GaussModelParam);
% Gauss1(amplParams(2:end)) = 0;
% Gauss2(amplParams([1 3])) = 0;
% Gauss3(amplParams(1:2)) = 0;
% G = GaussModel(GaussModelParam, freq(freqBounds));
% G1 = ThreeGaussModel_noBaseline(Gauss1, freq(freqBounds));
% G2 = ThreeGaussModel_noBaseline(Gauss2, freq(freqBounds));
% G3 = ThreeGaussModel_noBaseline(Gauss3, freq(freqBounds));

hold on;
plot(freq(freqBounds), L);
% plot(freq(freqBounds), L_baseline);
% plot(freq(freqBounds), G);
% plot(freq(freqBounds), G1);
% plot(freq(freqBounds), G2);
% plot(freq(freqBounds), G3);
% plot(freq(freqBounds), ThreeGaussModel_noBaseline(Gauss2, freq(freqBounds)));
% plot(freq(freqBounds), ThreeGaussModel_noBaseline(Gauss3, freq(freqBounds)));
legend({'data','model + baseline','baseline','residual','model'}, ...
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
exportgraphics(h_tmp, fullfile(out_dir, [fname '_Glu_model_fit.png']), "Resolution", 300);
close(h_tmp);

GluHeight = modelParam(1) * modelParam(4);
MRS_struct.out.(vox{kk}).Glu.FitError(ii)     = 100 * std(resid) / GluHeight;
MRS_struct.out.(vox{kk}).Glu.Area(ii)         = modelParam(1) * pi;
MRS_struct.out.(vox{kk}).Glu.FWHM(ii)         = 2 / modelParam(4) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).Glu.ModelParam(ii,:) = modelParam;
MRS_struct.out.(vox{kk}).Glu.Resid(ii,:)      = resid;

% Calculate SNR of Glu signal
noiseSigma_SUM = CalcNoise(freq, SUM(ii,:));
MRS_struct.out.(vox{kk}).Glu.SNR(ii) = abs(GluHeight) / noiseSigma_SUM;
