function [modelParam, rejectFrame, residCr] = FitChoCr(freq, spec, modelParamInit, LarmorFreq)

warning('off','stats:nlinfit:IterationLimitExceeded'); % temporarily suppress warning messages about iteration limit

% All parameters in initx are in standard units.
% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = [1 2*LarmorFreq LarmorFreq 180/pi 1 1 1];
modelParamInit = modelParamInit ./ conv;

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6);

n = size(spec,1);
modelParam = zeros(n,7);

for ii = 1:n
    modelParamInit = lsqcurvefit(@TwoLorentzModel, modelParamInit, freq, real(spec), [], [], lsqopts);
    [modelParam(ii,:), residCr] = nlinfit(freq, real(spec), @TwoLorentzModel, modelParamInit, nlinopts);
    
    % fit_plot = TwoLorentzModel(modelParam(ii,:), freq);
    % figure(3);
    % plot(freq, real(spec), 'g', freq, fit_plot,'b');
    % set(gca,'XDir','reverse');
    % drawnow;
end

for ii = 1:size(modelParam,1)
    if modelParam(ii,1) < 0
        modelParam(ii,4) = modelParam(ii,4) + pi;
    end
end

% Need to deal with phase wrap:
% Convert to complex number then recalculate phase within 2*pi range
phase_wrapped = modelParam(:,4);
phase_unwrapped = angle(complex(cos(phase_wrapped), sin(phase_wrapped)));

% then fix to be within -pi..pi
offsetpos =  pi*lt(phase_unwrapped, -pi/2);
offsetneg = -pi*gt(phase_unwrapped,  pi/2);
phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
modelParam(:,4) = phase_unwrapped;

% Fix linewidth to be positive
modelParam(:,2) = abs(modelParam(:,2));

% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = repmat([1 2*LarmorFreq LarmorFreq 180/pi 1 1 1], [n 1]);
modelParam = modelParam .* conv;

% Reject any point where the fit params - area, fwhm, phase
%  or freq are > 3stdev away from the mean
% Set reject criteria for all fit parameters
MeanFitParams = mean(modelParam, 1);
UpperLim = repmat(MeanFitParams + 3*std(modelParam,1), [n 1]);
LowerLim = repmat(MeanFitParams - 3*std(modelParam,1), [n 1]);
% But don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectFrame = gt(modelParam, UpperLim);
rejectFrame = rejectFrame + lt(modelParam, LowerLim);
rejectFrame = max(rejectFrame,[],2);

warning('on','stats:nlinfit:IterationLimitExceeded'); % turn warning about about iteration limit back on

end
