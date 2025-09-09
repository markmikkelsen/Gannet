function [FitParams, rejectframe, residCr]  = FitCr(freq, FrameData, initx, LarmorFreq)

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');

%All parameters in initx are in standard units.
% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = [1 2*LarmorFreq LarmorFreq 180/pi 1 1];
initx = initx./conv;

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',1e5,'Display','off');

nframes = size(FrameData,2);
FitParams = zeros(nframes,6);

for jj = 1:nframes
    initx = lsqcurvefit(@LorentzModel, initx, freq', real(FrameData(:,jj)), [], [], lsqopts);
    [FitParams(jj,:), residCr] = nlinfit(freq', real(FrameData(:,jj)), @LorentzModel, initx, nlinopts);
    
    %fit_plot = LorentzModel(fit_param, freq);
    %figure(3); plot(freq', real(FrameData(:,jj)), 'g', freq', fit_plot,'b');
    %pause(0.8)
    %set(gca,'XDir','reverse');
end

for kk = 1:size(FitParams,1)
    if FitParams(kk,1) < 0
        FitParams(kk,4) = FitParams(kk,4) + pi;
    end
end

% Need to deal with phase wrap:
% Convert to complex number then recalculate phase within 2*pi range
phase_wrapped = FitParams(:,4);
cmplx = complex(cos(phase_wrapped), sin(phase_wrapped));
phase_unwrapped = angle(cmplx);

% then fix to be within -pi..pi
offsetpos =  pi*lt(phase_unwrapped, -pi/2);
offsetneg = -pi*gt(phase_unwrapped,  pi/2);
phase_unwrapped = phase_unwrapped + offsetpos + offsetneg;
FitParams(:,4) = phase_unwrapped;

% Fix area and linewidth to be positive
FitParams(:,1) = abs(FitParams(:,1));
FitParams(:,2) = abs(FitParams(:,2));

% Conversion factors to FWHM in Hz, delta f0 in Hz, phase in degrees
conv = repmat([1 2*LarmorFreq LarmorFreq 180/pi 1 1], [nframes 1]);
FitParams = FitParams .* conv;

% Reject any point where the fit params - area, fwhm, phase
%  or freq are > 3stdev away from the mean
% set reject criteria for all fit parameters
MeanFitParams = mean(FitParams, 1);
UpperLim = repmat(MeanFitParams + 3*std(FitParams,1), [nframes 1]);
LowerLim = repmat(MeanFitParams - 3*std(FitParams,1), [nframes 1]);
%but don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectframe = gt(FitParams, UpperLim);
rejectframe = rejectframe + lt(FitParams, LowerLim);
rejectframe = max(rejectframe,[],2);

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','stats:nlinfit:IterationLimitExceeded');

end
