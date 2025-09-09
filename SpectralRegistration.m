function [AllFramesFTrealign, MRS_struct] = SpectralRegistration(MRS_struct, OnWhat, Dual)
% Spectral registration is a time-domain frequency-and-phase correction
% routine as per Near et al. (2015)
% OnWhat=0 for spectro data, OnWhat=1 for water data
% MM: updates to improve speed and robustness (Jun 2017)

ii = MRS_struct.ii;
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-8,'TolFun',1e-8);

% Dual-channel option only applies registration separately to ONs and OFFs
SpecRegLoop = 0;
if nargin == 3
    if Dual == 1
        % We want to run this code twice, once for ONs, once for OFFs.
        SpecRegLoop = 1;
    end
end

time = (0:1:(MRS_struct.p.npoints(ii)-1)).'/MRS_struct.p.sw(ii);

% Use first n points of time-domain data, where n is the last point where SNR > 1.5
signal = abs(MRS_struct.fids.data);
noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
SNR = signal ./ repmat(noise, [size(MRS_struct.fids.data,1) 1]);
SNR = mean(SNR,2);
n = find(SNR > 1.5);
if isempty(n)
    tMax = 100;
else
    tMax = n(end);
    if tMax < 100
        tMax = 100;
    end
end

while SpecRegLoop > -1
    
    if OnWhat %Read water data
        %First, take the complex data and turn it into a real matrix
        flatdata(:,1,:) = real(MRS_struct.fids.data_water(1:tMax,:));
        flatdata(:,2,:) = imag(MRS_struct.fids.data_water(1:tMax,:));
    else % read spectro data
        if nargin == 3
            if Dual == 1
                %This code runs twice, first for ONs, second for OFFs.
                %SpecRegLoop;
                %size(real(MRS_struct.fids.data(:,(MRS_struct.fids.ON_OFF==SpecRegLoop))));
                
                flatdata(:,1,:) = real(MRS_struct.fids.data(1:tMax,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
                flatdata(:,2,:) = imag(MRS_struct.fids.data(1:tMax,(MRS_struct.fids.ON_OFF==SpecRegLoop)));
            end
        else
            % First, take the complex data and turn it into a real matrix
            flatdata(:,1,:) = real(MRS_struct.fids.data(1:tMax,:));
            flatdata(:,2,:) = imag(MRS_struct.fids.data(1:tMax,:));
        end
    end
    
    %Correct to a point 10% into the file (seems better that the actual beginning)
    %AlignRow = ceil(size(flatdata,3)/10);
    %flattarget = squeeze(flatdata(:,:,AlignRow));
    
    % Use median across transients
    flattarget = median(flatdata,3); % median across transients
    
    %Time domain Frequency and Phase Correction
    %Preliminary to fitting:
    parsGuess = [0 0]; %initial freq and phase guess
    parsFit = zeros([size(flatdata,3) 2]);
    CorrPars = zeros([size(flatdata,3) 2]);
    MSE = zeros([1 size(flatdata,3)]);
    input.dwelltime = 1/MRS_struct.p.sw(ii);
    
    %Fitting to determine frequency and phase corrections
    target = flattarget(:);
    reverseStr = '';
    for corrloop = 1:size(flatdata,3)
        msg = sprintf('\nRunning spectral registration on transient: %d', corrloop);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        
        transient = squeeze(flatdata(:,:,corrloop));
        input.data = transient(:);
        [parsFit(corrloop,:), ~, ~, ~, MSE(corrloop)] = nlinfit(input, target, @FreqPhaseShiftNest, parsGuess, nlinopts);
        parsGuess = parsFit(corrloop,:); %Carry parameters from point to point
    end
    zMSE = zscore(MSE); % standardized MSEs
    
    %Applying frequency and phase corrections
    if OnWhat        
        for corrloop=1:size(flatdata,3)
            MRS_struct.fids.data_water(:,corrloop)=MRS_struct.fids.data_water(:,corrloop).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
        end
        FullData = MRS_struct.fids.data_water;
        FullData = FullData.* repmat( (exp(-(time)*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data_water,2)]);
        AllFramesFTrealign = fftshift(fft(FullData, MRS_struct.p.ZeroFillTo(ii),1),1);
    else
        MRS_struct.out.SpecReg.freq{ii}  = parsFit(:,1);
        MRS_struct.out.SpecReg.phase{ii} = parsFit(:,2);
        
        for corrloop=1:size(flatdata,3)
            if nargin == 3
                if Dual == 1
                    %Need to get the slot right to put data back into
                    averages_per_dynamic=find(MRS_struct.fids.ON_OFF~=(MRS_struct.fids.ON_OFF(1)),1)-1;
                    dyn=floor((corrloop-1)/averages_per_dynamic); %number of cycles in
                    ind=mod((corrloop-1),averages_per_dynamic)+1; %number in current cycle
                    
                    if SpecRegLoop == 1
                        if MRS_struct.fids.ON_OFF(1) == 1
                            corrloop_d = dyn*averages_per_dynamic*2+ind;
                        else
                            corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                        end
                    else
                        if MRS_struct.fids.ON_OFF(1) == 1
                            corrloop_d = dyn*averages_per_dynamic*2+averages_per_dynamic+ind;
                        else
                            corrloop_d = dyn*averages_per_dynamic*2+ind;
                        end
                    end
                    
                    MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
                end
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            else
                corrloop_d=corrloop;
                MRS_struct.fids.data_align(:,corrloop_d)=MRS_struct.fids.data(:,corrloop_d).*exp(1i*parsFit(corrloop,1)*2*pi*time)*exp(1i*pi/180*parsFit(corrloop,2));
                CorrPars(corrloop_d,:)=parsFit(corrloop,:);
            end
        end
        
        if SpecRegLoop == 0
            FullData = MRS_struct.fids.data_align;
            FullData = FullData .* repmat((exp(-(time) * MRS_struct.p.LB * pi)), [1 size(MRS_struct.fids.data,2)]);
            AllFramesFTrealign = fftshift(fft(FullData, MRS_struct.p.ZeroFillTo(ii),1),1);
            
            % In frequency domain, move Cr to 3.02 and get phase 'right' as opposed to 'consistent'
            if ~MRS_struct.p.phantom
                freqbounds = MRS_struct.spec.freq <= 3.6 & MRS_struct.spec.freq >= 2.6; % MM (170227)
                freq = MRS_struct.spec.freq(freqbounds);
                %Do some detective work to figure out the initial parameters
                ChoCrMeanSpec = mean(AllFramesFTrealign(freqbounds,:),2);
                Baseline_offset = real(ChoCrMeanSpec(1)+ChoCrMeanSpec(end))/2;
                Width_estimate = 0.05; %ppm
                Area_estimate = (max(real(ChoCrMeanSpec))-min(real(ChoCrMeanSpec)))*Width_estimate*4;
                ChoCr_initx = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
                
                if nargin == 3
                    if Dual == 1
                        %This bit is silly - we don't want to do OFF-to-ON based on the Cr signal
                        ChoCrMeanSpecON = mean(AllFramesFTrealign(freqbounds,(MRS_struct.fids.ON_OFF==1)),2);
                        ChoCrMeanSpecOFF = mean(AllFramesFTrealign(freqbounds,(MRS_struct.fids.ON_OFF==0)),2);
                        ChoCrMeanSpecFitON = FitChoCr(freq, ChoCrMeanSpecON, ChoCr_initx,MRS_struct.p.LarmorFreq(ii));
                        ChoCrMeanSpecFitOFF = FitChoCr(freq, ChoCrMeanSpecOFF, ChoCr_initx,MRS_struct.p.LarmorFreq(ii));
                        AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1))*exp(1i*pi/180*(ChoCrMeanSpecFitON(4))); % phase
                        AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))=AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0))*exp(1i*pi/180*(ChoCrMeanSpecFitOFF(4))); % phase
                        
                        ChoCrFreqShiftON = ChoCrMeanSpecFitON(3);
                        ChoCrFreqShiftON = ChoCrFreqShiftON - 3.02*MRS_struct.p.LarmorFreq(ii);
                        ChoCrFreqShiftON = ChoCrFreqShiftON ./ (MRS_struct.p.LarmorFreq(ii) * abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
                        ChoCrFreqShift_pointsON = round(ChoCrFreqShiftON);
                        AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)) = circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==1)), ChoCrFreqShift_pointsON); % freq
                        ChoCrFreqShiftOFF = ChoCrMeanSpecFitOFF(3);
                        ChoCrFreqShiftOFF = ChoCrFreqShiftOFF - 3.02*MRS_struct.p.LarmorFreq(ii);
                        ChoCrFreqShiftOFF = ChoCrFreqShiftOFF ./ (MRS_struct.p.LarmorFreq(ii) * abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
                        ChoCrFreqShift_pointsOFF = round(ChoCrFreqShiftOFF);
                        AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)) = circshift(AllFramesFTrealign(:,(MRS_struct.fids.ON_OFF==0)), ChoCrFreqShift_pointsOFF); % freq
                        
                    end
                else
                    ChoCrMeanSpecFit = FitChoCr(freq, ChoCrMeanSpec, ChoCr_initx, MRS_struct.p.LarmorFreq(ii));
                    AllFramesFTrealign = AllFramesFTrealign*exp(1i*pi/180*(ChoCrMeanSpecFit(4))); % phase
                    ChoCrFreqShift = ChoCrMeanSpecFit(3);
                    ChoCrFreqShift = ChoCrFreqShift - 3.02*MRS_struct.p.LarmorFreq(ii);
                    ChoCrFreqShift = ChoCrFreqShift ./ (MRS_struct.p.LarmorFreq(ii) * abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
                    ChoCrFreqShift_points = round(ChoCrFreqShift);
                    AllFramesFTrealign = circshift(AllFramesFTrealign, ChoCrFreqShift_points); % freq
                end
            end
            
            %Some output
            MRS_struct.out.FreqStdevHz(ii) = std(parsFit(:,1),1);
            
            % Reject transients that are greater than 3 st. devs. of zMSEs
            % (only applies if not using weighted averaging)
            MRS_struct.out.SpecReg.zMSE{ii} = zMSE;
            MRS_struct.out.reject{ii}       = zMSE > 3;
            
        end        
    end
    
    SpecRegLoop = SpecRegLoop-1;
    
end

fprintf('\n');

end



