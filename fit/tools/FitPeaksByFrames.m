function [modelParam, rejectFrame, resid] = FitPeaksByFrames(freq, spec, modelParamInit)

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6,'FunValCheck','off');

n = size(spec,1);
modelParam = zeros(n,6);

for ii = 1:n
    modelParamInit = lsqcurvefit(@LorentzModel, modelParamInit, freq, real(spec(ii,:)), [], [], lsqopts);
    [modelParam(ii,:), resid] = nlinfit(freq, real(spec(ii,:)), @LorentzModel, modelParamInit, nlinopts);
        
    % fit_plot = LorentzModel(FitParams(ii,:), freq);    
    % figure(3);
    % plot(freq, real(FrameData(:,ii)), 'g', freq, fit_plot, 'b');
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
offsetpos =  2*pi*lt(phase_unwrapped, -pi);
offsetneg = -2*pi*gt(phase_unwrapped,  pi);
phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
modelParam(:,4) = phase_unwrapped;

% Fix area and linewidth to be positive
modelParam(:,1) = abs(modelParam(:,1));
modelParam(:,2) = abs(modelParam(:,2));

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
