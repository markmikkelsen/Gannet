function [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign, MRS_struct)

switch MRS_struct.p.alignment
    case 'NAA'
        peak_ppm = 2.02;
        ppm_lim = MRS_struct.spec.freq > 1.76 & MRS_struct.spec.freq < 2.26;
    case 'Cr'
        peak_ppm = 3.02;
        ppm_lim = MRS_struct.spec.freq > 2.72 & MRS_struct.spec.freq < 3.12;
    case 'Cho'
        peak_ppm = 3.2;
        ppm_lim = MRS_struct.spec.freq > 3.1 & MRS_struct.spec.freq < 3.32;
end
freqRange = MRS_struct.spec.freq(ppm_lim);

% Set initial parameters by fitting the SUM
SUM = mean(AllFramesFTrealign(ppm_lim,:),2);
area = max(real(SUM)) / (2*pi);
modelParamInit = [area 0.03 peak_ppm 0 0 0];
modelParam = FitPeaksByFrames(freqRange, SUM.', modelParamInit);

% Update phase and freq of AllFramesFTrealign to reflect SUM fit
AllFramesFTrealign = AllFramesFTrealign * exp(1i*modelParam(4));
freq_step_size = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
shift_points = round((modelParam(3) - peak_ppm) / freq_step_size);
for ii = 1:size(AllFramesFTrealign,2)
    AllFramesFTrealign(:,ii) = circshift(AllFramesFTrealign(:,ii), [shift_points 0]);
end

% Set initial fitting parameters by fitting each average
modelParamInit = modelParam;
modelParamInit(3:4) = [peak_ppm 0];
modelParam = FitPeaksByFrames(freqRange, AllFramesFTrealign(ppm_lim,:).', modelParamInit);

% Update phase and freq of AllFramesFTrealign to reflect frame-by-frame fit
AllFramesFTrealign = AllFramesFTrealign .* repmat(exp(1i*modelParam(:,4)).', [length(MRS_struct.spec.freq) 1]);
shift_points = round((modelParam(:,3) - peak_ppm) / freq_step_size);
for ii = 1:size(shift_points,1)
    AllFramesFTrealign(:,ii) = circshift(AllFramesFTrealign(:,ii),[shift_points(ii) 0]);
end

% Fit just the Cr in the aligned mean spectrum to get CrFWHMHz
ppm_lim = MRS_struct.spec.freq > 2.6 & MRS_struct.spec.freq < 3.11;
freqRange = MRS_struct.spec.freq(ppm_lim);
CrMeanSpec = mean(AllFramesFTrealign(ppm_lim,:),2);
CrMeanSpecFit = FitCr(freqRange, CrMeanSpec.', modelParamInit, MRS_struct.p.LarmorFreq(MRS_struct.ii));
MRS_struct.out.CrFWHMHz(MRS_struct.ii) = CrMeanSpecFit(2);

% Reject any point where the fit params - area, FWHM, phase
%  or freq are >3 stdev away from the mean
% Set reject criteria for all fit parameters
MeanFrameParams = mean(modelParam, 1);
UpperLim = repmat(MeanFrameParams + 3*std(modelParam,1), [size(AllFramesFTrealign,2) 1]);
LowerLim = repmat(MeanFrameParams - 3*std(modelParam,1), [size(AllFramesFTrealign,2) 1]);
% But don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectframe = gt(modelParam, UpperLim);
rejectframe = rejectframe + lt(modelParam, LowerLim);
MRS_struct.out.reject{MRS_struct.ii} = max(rejectframe,[],2)';

% Balance up rejects
stepsize = find(MRS_struct.fids.ON_OFF ~= MRS_struct.fids.ON_OFF(1), 1) - 1;
for ii = 1:length(MRS_struct.out.reject{MRS_struct.ii})
    % First find whether reject is ON or OFF
    if MRS_struct.out.reject{MRS_struct.ii}(ii) == 1
        if MRS_struct.fids.ON_OFF(ii) == MRS_struct.fids.ON_OFF(1)
            MRS_struct.out.reject{MRS_struct.ii}(ii + stepsize) = 1;
        else
            MRS_struct.out.reject{MRS_struct.ii}(ii - stepsize) = 1;
        end
    end
end

fprintf('\n');

end
