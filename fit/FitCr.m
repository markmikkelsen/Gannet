function [modelParam, rejectFrame, residCr] = FitCr(freq, spec, modelParamInit, LarmorFreq)

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');

% All parameters in initx are in standard units.
% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = [1 2*LarmorFreq LarmorFreq 180/pi 1 1];
modelParamInit = modelParamInit ./ conv;

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',1e5,'Display','off');

n = size(spec,1);
modelParam = zeros(n,6);

for ii = 1:n
    modelParamInit = lsqcurvefit(@LorentzModel, modelParamInit, freq, real(spec(ii,:)), [], [], lsqopts);
    [modelParam(ii,:), residCr] = nlinfit(freq, real(spec(ii,:)), @LorentzModel, modelParamInit, nlinopts);
    
    % fit_plot = LorentzModel(modelParam, freq);
    % figure(3); plot(freq, real(spec(ii,:)), 'g', freq, fit_plot, 'b');
    % pause(0.8);
    % set(gca,'XDir','reverse');
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

% Fix area and linewidth to be positive
modelParam(:,1) = abs(modelParam(:,1));
modelParam(:,2) = abs(modelParam(:,2));

% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = repmat([1 2*LarmorFreq LarmorFreq 180/pi 1 1], [n 1]);
modelParam = modelParam .* conv;

% Reject any point where the fit params - area, fwhm, phase
%  or freq are > 3stdev away from the mean
% Set reject criteria for all fit parameters
MeanFitParams = mean(modelParam,1);
UpperLim = repmat(MeanFitParams + 3*std(modelParam,1), [n 1]);
LowerLim = repmat(MeanFitParams - 3*std(modelParam,1), [n 1]);
% But don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectFrame = gt(modelParam, UpperLim);
rejectFrame = rejectFrame + lt(modelParam, LowerLim);
rejectFrame = max(rejectFrame,[],2);

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','stats:nlinfit:IterationLimitExceeded');

end
