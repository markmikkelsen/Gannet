function [MRS_struct, Cr_OFF, freqBoundsChoCr, freqBoundsCr, residCr] = FitCrNAAGlu(MRS_struct, freq, OFF, SUM, vox, ii, kk, lsqopts, nlinopts)
                
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
LorentzModelInit = [max(real(Cr_OFF(freqBoundsCr))) 0.05 3.0 0 0 0];

% Least-squares model fitting
LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqBoundsCr), real(Cr_OFF(freqBoundsCr)), [], [], lsqopts);
[LorentzModelParam, residCr] = nlinfit(freq(freqBoundsCr), real(Cr_OFF(freqBoundsCr)), @LorentzModel, LorentzModelInit, nlinopts);

MRS_struct.out.(vox{kk}).Cr.ModelParam(ii,:) = LorentzModelParam;
CrHeight = LorentzModelParam(1) / (2*pi*LorentzModelParam(2));
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

LorentzModelInit = [maxinNAA 0.05 2.01 0 -LinearInit constInit];
lb = [0 0.01 1.97 0 -40*maxinNAA -2000*maxinNAA];
ub = [4000*maxinNAA 0.1 2.05 0.5 40*maxinNAA 1000*maxinNAA];

% Least-squares model fitting
LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqBounds), real(NAA_OFF(freqBounds)), lb, ub, lsqopts);
[LorentzModelParam, resid] = nlinfit(freq(freqBounds), real(NAA_OFF(freqBounds)), @LorentzModel, LorentzModelInit, nlinopts);

NAAheight = LorentzModelParam(1) / (2*pi*LorentzModelParam(2));
MRS_struct.out.(vox{kk}).NAA.FitError(ii) = 100 * std(resid) / NAAheight;
NAAModelParam = LorentzModelParam;
NAAModelParam(4) = 0;
MRS_struct.out.(vox{kk}).NAA.Area(ii) = sum(LorentzModel(NAAModelParam,freq(freqBounds)) - BaselineModel(NAAModelParam([3 6 5]),freq(freqBounds)), 2) * abs(freq(1) - freq(2));
MRS_struct.out.(vox{kk}).NAA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * NAAModelParam(2));
MRS_struct.out.(vox{kk}).NAA.ModelParam(ii,:) = LorentzModelParam;
MRS_struct.out.(vox{kk}).NAA.Resid(ii,:) = resid;

% Calculate SNR of NAA signal
MRS_struct.out.(vox{kk}).NAA.SNR(ii) = abs(NAAheight) / noiseSigma_OFF;


% Glu Fit - MM (240328): This still needs testing to see if it's reliable

Glu_SUM    = SUM(ii,:);
freqBounds = find(freq <= 2.43 & freq >= 2.25);

maxinGlu    = max(real(Glu_SUM(freqBounds)));
gradPoints = (real(Glu_SUM(freqBounds(end))) - real(Glu_SUM(freqBounds(1)))) ./ abs(freqBounds(end) - freqBounds(1)); %in points
LinearInit  = gradPoints ./ abs(freq(1) - freq(2));
constInit   = (real(Glu_SUM(freqBounds(end))) + real(Glu_SUM(freqBounds(1)))) ./ 2;

LorentzModelInit = [maxinGlu -0.2*maxinGlu -0.2*maxinGlu ...
                    30 ...
                    2.34 ...
                    0.05 ...
                    0 ...
                   -LinearInit -LinearInit constInit];

lb = [-4e3*maxinGlu -800*maxinGlu -800*maxinGlu ...
       0 ...
       2.34-0.05 ...
       0 ...
      -180 ...
      -40*maxinGlu -40*maxinGlu -2000*maxinGlu];

ub = [4e3*maxinGlu 800*maxinGlu 800*maxinGlu ...
      100 ...
      2.34+0.05 ...
      0.1 ...
      180 ...
      40*maxinGlu 40*maxinGlu 1000*maxinGlu];

% Least-squares model fitting
LorentzModelInit = lsqcurvefit(@ThreeLorentzModel, LorentzModelInit, freq(freqBounds), real(Glu_SUM(freqBounds)), lb, ub, lsqopts);
[LorentzModelParam, resid] = nlinfit(freq(freqBounds), real(Glu_SUM(freqBounds)), @ThreeLorentzModel, LorentzModelInit, nlinopts);

GluHeight = LorentzModelParam(1) * LorentzModelParam(4);
MRS_struct.out.(vox{kk}).Glu.FitError(ii)     = 100 * std(resid) / GluHeight;
MRS_struct.out.(vox{kk}).Glu.Area(ii)         = LorentzModelParam(1) * pi;
MRS_struct.out.(vox{kk}).Glu.FWHM(ii)         = 2 / LorentzModelParam(4) * MRS_struct.p.LarmorFreq(ii);
MRS_struct.out.(vox{kk}).Glu.ModelParam(ii,:) = LorentzModelParam;
MRS_struct.out.(vox{kk}).Glu.Resid(ii,:)      = resid;

% Calculate SNR of Glu signal
noiseSigma_SUM = CalcNoise(freq, SUM(ii,:));
MRS_struct.out.(vox{kk}).Glu.SNR(ii) = abs(GluHeight) / noiseSigma_SUM;
