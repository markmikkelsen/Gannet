function [AllFramesFTrealign, MRS_struct] = AlignUsingPeak(AllFramesFTrealign, MRS_struct)
%Using the freq information from MRS_struct, determine appropriate range to
%fit NAA signal over
switch MRS_struct.p.alignment
    case 'NAA'
        PEAK_ppm = 2.02;
        %Determine limits
        z = abs(MRS_struct.spec.freq-2.26);
        lb = find(min(z) == z);
        z = abs(MRS_struct.spec.freq-1.76);
        ub = find(min(z) == z);
        Initx = [50 0.1 PEAK_ppm 0 0 0];
        
    case 'Cr'
        PEAK_ppm = 3.02;
        %Determine limits
        z = abs(MRS_struct.spec.freq-3.12);
        lb = find(min(z) == z);
        z = abs(MRS_struct.spec.freq-2.72);
        ub = find(min(z) == z);
        Initx = [30 0.05 PEAK_ppm 0 0 0];
        
    case 'Cho'
        PEAK_ppm = 3.2;
        %Determine limits
        z = abs(MRS_struct.spec.freq-3.32);
        lb = find(min(z) == z);
        z = abs(MRS_struct.spec.freq-3.1);
        ub = find(min(z) == z);
        Initx = [30 0.05 PEAK_ppm 0 0 0];
end
%Set initial parameters by fitting the sum

freqrange = MRS_struct.spec.freq(lb:ub);
SumSpec = sum(AllFramesFTrealign(lb:ub,:),2);
SumSpecParams = FitPeaksByFrames2(freqrange, SumSpec, Initx);
%FitPeaksByFrames2 now edited so that output params have same units as Initx
%Update freq and phase of AllFramesFTrealign to reflect Sum fit
AllFramesFTrealign = AllFramesFTrealign * exp(1i*SumSpecParams(4));
%Shift NAA freq to 2.02
freq_step_size = abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2));
Shift_points = round((real(SumSpecParams(3)) - PEAK_ppm)/freq_step_size);
for ii = 1:size(AllFramesFTrealign,2)
    AllFramesFTrealign(:,ii) = circshift(AllFramesFTrealign(:,ii), [Shift_points 0]);
end
%Set initial fitting parameters by fitting the sum
Init = SumSpecParams ./ [size(AllFramesFTrealign,2) 1 1 1 size(AllFramesFTrealign,2) size(AllFramesFTrealign,2)];
Init(3) = PEAK_ppm;
Init(4) = 0;
%Select Data2BeFit
%figure(7)
%subplot(1,3,1)
%plot(MRS_struct.spec.freq,real(AllFramesFTrealign))
%set(gca,'XLim',[1.5 2.5]);
%set(gca,'XDir','reverse')
Data2bFit = AllFramesFTrealign(lb:ub,:);
FrameParams = FitPeaksByFrames2(freqrange, Data2bFit, Init); %FitPeaksByFrames2 now edited so that output params have same units as Initx
%Update freq and phase of AllFramesFTrealign to reflect Frame-by-Frame fit
FrameShift_points = round((real(FrameParams(:,3))-PEAK_ppm)/freq_step_size);
for jj = 1:size(FrameShift_points,1)
    AllFramesFTrealign(:,jj) = circshift(AllFramesFTrealign(:,jj),[FrameShift_points(jj) 0]);
end
%subplot(1,3,2)
%plot(MRS_struct.spec.freq,real(AllFramesFTrealign))
%set(gca,'XLim',[1.5 2.5]);
%set(gca,'XDir','reverse')
%mean(FrameShift_points)
%mean(FrameParams(:,3))
%size(AllFramesFTrealign)
%size(repmat(exp(1i*SumSpecParams(:,4)).',[length(MRS_struct.spec.freq) 1]))
AllFramesFTrealign = AllFramesFTrealign.*repmat(exp(1i*FrameParams(:,4)).',[length(MRS_struct.spec.freq) 1]);
AllFramesFTrealign = AllFramesFTrealign-repmat(FrameParams(:,5).',[length(MRS_struct.spec.freq) 1]);
%subplot(1,3,3)
%plot(MRS_struct.spec.freq,real(AllFramesFTrealign))
%set(gca,'XLim',[1.5 2.5]);
%set(gca,'XDir','reverse');
%Fit just the Cr in the aligned mean spectrum to get CrFWHMHz

CrFitLimLow = 2.6;
CrFitLimHigh = 3.11;
%Still need ranges for Creatine align plot
z = abs(MRS_struct.spec.freq-CrFitLimHigh);
clb = find(min(z) == z);
z = abs(MRS_struct.spec.freq-CrFitLimLow);
cub = find(min(z) == z);
freqrange = MRS_struct.spec.freq(clb:cub);
%Cr_initx = [ Area_estimate Width_estimate 3.02 0 Baseline_offset 0 ].*[1 (2*MRS_struct.p.LarmorFreq) (MRS_struct.p.LarmorFreq) (180/pi) 1 1 ];
CrMeanSpec = mean(AllFramesFTrealign(clb:cub,:),2);
CrMeanSpecFit = FitCr(freqrange, CrMeanSpec, Init, MRS_struct.p.LarmorFreq);

%Some Output
MRS_struct.out.CrFWHMHz(MRS_struct.ii) = CrMeanSpecFit(2);

%plot the models.
%     Models=zeros(size(Data2bFit));
%     for ii=1:size(Models,2)
%         Models(:,ii)=LorentzModel(FrameParams(ii,:),freqrange);
%     end
%     subplot(3,1,1)
%     plot(freqrange,Models);
%     set(gca,'XLim',[1.5 2.5]);
%     set(gca,'XDir','reverse');
%     for ii=1:size(Models,2)
%     %Models_undone(:,ii)=((Models(:,ii)-FrameParams(ii,5)-FrameParams(ii,6)*(freqrange.'-FrameParams(ii,3))))*exp(-1i*FrameParams(ii,4));
%     Models_undone(:,ii)=(Models(:,ii)*exp(1i*FrameParams(ii,4)));
%     %Models_undone(:,ii)=circshift(Models_undone(:,ii),[FrameShift_points(jj) 0]);
%     end
%     subplot(3,1,2)
%     plot(freqrange,Models_undone);
%     set(gca,'XLim',[1.5 2.5]);
%     set(gca,'XDir','reverse');
%    subplot(3,1,3)
%    plot(1:length(freqrange),std(Models.'),1:length(freqrange),std(Models_undone.'))
%    plot(freqrange,sum(Data2bFit,2),freqrange,sum(Models_undone,2))
%     save testing.mat Data2bFit FrameParams freqrange
%plot(FrameParams(:,4)/pi*180)


% Reject any point where the fit params - area, fwhm, phase
%  or freq are > 3stdev away from the mean
% set reject criteria for all fit parameters
MeanFrameParams = mean(FrameParams, 1);
UpperLim = repmat(MeanFrameParams + 3*std(FrameParams,1), [size(AllFramesFTrealign,2) 1]);
LowerLim = repmat(MeanFrameParams - 3*std(FrameParams,1), [size(AllFramesFTrealign,2) 1]);
%but don't reject on linear, const baseline fit vals
UpperLim(:,5:6) = Inf;
LowerLim(:,5:6) = -Inf;
rejectframe = gt(FrameParams, UpperLim);
rejectframe = rejectframe + lt(FrameParams, LowerLim);
MRS_struct.out.reject(:,MRS_struct.ii) = max(rejectframe,[],2);
%Balance up rejects

IsRejectOnOrOff = MRS_struct.fids.ON_OFF(MRS_struct.out.reject(:,MRS_struct.ii) == 1);
stepsize = find(MRS_struct.fids.ON_OFF ~= (MRS_struct.fids.ON_OFF(1)),1)-1;
for kk = 1:size(MRS_struct.out.reject(:,MRS_struct.ii),1)
    %first find whether reject is ON or OFF
    if(MRS_struct.out.reject(kk,MRS_struct.ii) == 1)
        IsRejectOnOrOff = MRS_struct.fids.ON_OFF(kk);
        if IsRejectOnOrOff == MRS_struct.fids.ON_OFF(1)
            MRS_struct.out.reject(kk+stepsize, MRS_struct.ii) = 1;
        else
            MRS_struct.out.reject(kk-stepsize, MRS_struct.ii) = 1;
        end
    end
end

end



