function [AllFramesFTrealign, MRS_struct] = SpectralRegistrationHERMES(MRS_struct)
% Align using a multistep variant of spectral registration (Mikkelsen et
% al. Magn Reson Med. 2018;80(1):21-28. doi:10.1002/mrm.27027)

showPlots = 0;

% Looping parameters
if MRS_struct.p.HERMES % run registration four times - once for each HERMES experiment
    SpecRegLoop = 3;
    SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
else % run registration once or twice for MEGA-PRESS acquisitions
    SpecRegLoop = 1;
    SubspecToAlign = MRS_struct.fids.ON_OFF;
end

% Pre-allocate memory
ii = MRS_struct.ii;
MRS_struct.out.SpecReg.freq{ii}  = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.phase{ii} = zeros(1,size(MRS_struct.fids.data,2));
zMSE = zeros(1,size(MRS_struct.fids.data,2));
CorrParsML = zeros(size(MRS_struct.fids.data,2),2);
count = 0;
parsGuess = [0 0];

% Inputs
DataToAlign = MRS_struct.fids.data;
time = (0:1:(MRS_struct.p.npoints(ii)-1)).'/MRS_struct.p.sw(ii);
input.dwelltime = 1/MRS_struct.p.sw(ii);

% Probability density function and parameter bounds
Cauchy = @(x,s,l) s./(pi.*(s.^2+(x-l).^2));
lb = [0 -Inf];
ub = [Inf Inf];

% Optimization options
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-8,'TolFun',1e-8);
mleopts  = statset('mlecustom');
mleopts  = statset(mleopts,'MaxIter',400,'MaxFunEvals',800,'TolX',1e-6,'TolFun',1e-6,'TolBnd',1e-6);

% Set dimensions of figures of histograms
if showPlots == 1
    d.w = 0.6;
    d.h = 0.45;
    d.l = (1-d.w)/2;
    d.b = (1-d.h)/2;
end

count2 = 1;
reverseStr = '';
while SpecRegLoop > -1
    
    % Use first n points of time-domain data, where n is the last point where SNR > 1.5
    signal = abs(DataToAlign(:,SubspecToAlign == SpecRegLoop));
    noise = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
    SNR = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
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
    
    % 'Flatten' complex data for use in nlinfit
    clear flatdata;
    flatdata(:,1,:) = real(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    
    % Reference transient
    flattarget = median(flatdata,3); % median across transients
    target = flattarget(:);
    
    % Pre-allocate memory
    if ~count
        parsFit = zeros(size(flatdata,3), 2);
        MSE = zeros(1, size(flatdata,3));
    end
    
    % Determine frequency and phase offsets by spectral registration
    for corrloop = 1:size(flatdata,3)
        msg = sprintf('\nRunning spectral registration (HERMES) on transient: %d', count2);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        count2 = count2 + 1;
        
        transient = squeeze(flatdata(:,:,corrloop));
        input.data = transient(:);
        [parsFit(corrloop,:), ~, ~, ~, MSE(corrloop)] = nlinfit(input, target, @FreqPhaseShiftNest, parsGuess, nlinopts);
        parsGuess = parsFit(corrloop,:);
    end
    
    count = count + 1;
    
    % Probability distribution of frequency offsets (estimated by maximum likelihood)
    MRS_struct.out.MLalign.f.x(count,:,ii) = parsFit(:,1);
    start = [iqr(MRS_struct.out.MLalign.f.x(count,:,ii))/2, median(MRS_struct.out.MLalign.f.x(count,:,ii))];
    [MRS_struct.out.MLalign.f.p(count,:,ii), MRS_struct.out.MLalign.f.p_ci(:,:,count,ii)] = ...
        mle(MRS_struct.out.MLalign.f.x(count,:,ii), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    MRS_struct.out.MLalign.f.fx(count,:,ii) = ...
        linspace(1.5*min(MRS_struct.out.MLalign.f.x(count,:,ii)), 1.5*max(MRS_struct.out.MLalign.f.x(count,:,ii)), 1e3);
    MRS_struct.out.MLalign.f.pdf(count,:,ii) = Cauchy(MRS_struct.out.MLalign.f.fx(count,:,ii), ...
        MRS_struct.out.MLalign.f.p(count,1,ii), MRS_struct.out.MLalign.f.p(count,2,ii));
    
    % Probability distribution of phase offsets (estimated by maximum likelihood)
    MRS_struct.out.MLalign.ph.x(count,:,ii) = parsFit(:,2);
    start = [iqr(MRS_struct.out.MLalign.ph.x(count,:,ii))/2, median(MRS_struct.out.MLalign.ph.x(count,:,ii))];
    [MRS_struct.out.MLalign.ph.p(count,:,ii), MRS_struct.out.MLalign.ph.p_ci(:,:,count,ii)] = ...
        mle(MRS_struct.out.MLalign.ph.x(count,:,ii), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
    MRS_struct.out.MLalign.ph.fx(count,:,ii) = ...
        linspace(1.5*min(MRS_struct.out.MLalign.ph.x(count,:,ii)), 1.5*max(MRS_struct.out.MLalign.ph.x(count,:,ii)), 1e3);
    MRS_struct.out.MLalign.ph.pdf(count,:,ii) = Cauchy(MRS_struct.out.MLalign.ph.fx(count,:,ii), ...
        MRS_struct.out.MLalign.ph.p(count,1,ii), MRS_struct.out.MLalign.ph.p(count,2,ii));
    
    if showPlots == 1
        % Histogram of frequency offsets
        H1 = figure(333);
        set(H1, 'Color', 'w', 'Units', 'Normalized', 'OuterPosition', [d.l d.b d.w d.h]);
        subplot(1,2,1);
        bins = linspace(min(MRS_struct.out.MLalign.f.x(count,:,ii)), max(MRS_struct.out.MLalign.f.x(count,:,ii)), 15);
        binWidth = abs(bins(1)-bins(2));
        h = bar(bins, histcounts(MRS_struct.out.MLalign.f.x(count,:,ii), length(bins))/(length(MRS_struct.out.MLalign.f.x(count,:,ii))*binWidth), 'histc');
        h.FaceColor = [0.8 0.8 0.8];
        hold on;
        plot(MRS_struct.out.MLalign.f.fx(count,:,ii), MRS_struct.out.MLalign.f.pdf(count,:,ii), 'Color', [1 0 0], 'LineWidth', 1.2);
        hold off;
        xlabel('\Deltaf (Hz)', 'FontSize', 15);
        ylabel('P(x)', 'FontSize', 15);
        set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
        
        % Histogram of phase offsets
        subplot(1,2,2);
        bins = linspace(min(MRS_struct.out.MLalign.ph.x(count,:,ii)), max(MRS_struct.out.MLalign.ph.x(count,:,ii)), 15);
        binWidth = abs(bins(1)-bins(2));
        h = bar(bins, histcounts(MRS_struct.out.MLalign.ph.x(count,:,ii), length(bins))/(length(MRS_struct.out.MLalign.ph.x(count,:,ii))*binWidth), 'histc');
        h.FaceColor = [0.8 0.8 0.8];
        hold on
        plot(MRS_struct.out.MLalign.ph.fx(count,:,ii), MRS_struct.out.MLalign.ph.pdf(count,:,ii), 'Color', [1 0 0], 'LineWidth', 1.2)
        hold off
        xlabel('\Delta\phi (deg)', 'FontSize', 15);
        ylabel('P(x)', 'FontSize', 15);
        set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
        
        drawnow;
        %pause(1);
    end
    
    corrloop_d = find(SubspecToAlign == SpecRegLoop);
    MRS_struct.out.SpecReg.freq{ii}(corrloop_d)  = parsFit(:,1) - MRS_struct.out.MLalign.f.p(count,2,ii)';
    MRS_struct.out.SpecReg.phase{ii}(corrloop_d) = parsFit(:,2) - MRS_struct.out.MLalign.ph.p(count,2,ii)';
    CorrParsML(corrloop_d,1) = parsFit(:,1) - MRS_struct.out.MLalign.f.p(count,2,ii)';
    CorrParsML(corrloop_d,2) = parsFit(:,2) - MRS_struct.out.MLalign.ph.p(count,2,ii)';
    zMSE(corrloop_d) = zscore(MSE); % standardized MSEs
    
    % Apply frequency and phase corrections
    for corrloop = 1:size(flatdata,3)
        % Default correction
        %DataToAlign(:,corrloop_d(corrloop)) = DataToAlign(:,corrloop_d(corrloop)) .* ...
        %    exp(1i*parsFit(corrloop,1)*2*pi*time) * exp(1i*pi/180*parsFit(corrloop,2));
        
        % Freq/phase correction + Cauchy pdf location parameter shift
        DataToAlign(:,corrloop_d(corrloop)) = DataToAlign(:,corrloop_d(corrloop)) .* ...
            exp(1i*(parsFit(corrloop,1) - MRS_struct.out.MLalign.f.p(count,2,ii))*2*pi*time) * ...
            exp(1i*pi/180*(parsFit(corrloop,2) - MRS_struct.out.MLalign.ph.p(count,2,ii)));
    end
    
    if SpecRegLoop == 0
        
        if showPlots == 1            
            MRS_struct.out.MLalign.f_aligned.x(ii,:) = CorrParsML(:,1);
            start = [std(MRS_struct.out.MLalign.f_aligned.x(ii,:))/2, median(MRS_struct.out.MLalign.f_aligned.x(ii,:))];
            MRS_struct.out.MLalign.f_aligned.p(ii,:) = ...
                mle(MRS_struct.out.MLalign.f_aligned.x(ii,:), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
            MRS_struct.out.MLalign.f_aligned.fx(ii,:) = ...
                linspace(1.1*min(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 1.1*max(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 1e3);
            MRS_struct.out.MLalign.f_aligned.pdf(ii,:) = ...
                Cauchy(MRS_struct.out.MLalign.f_aligned.fx(ii,:), MRS_struct.out.MLalign.f_aligned.p(ii,1), MRS_struct.out.MLalign.f_aligned.p(ii,2));
            
            MRS_struct.out.MLalign.ph_aligned.x(ii,:) = CorrParsML(:,2);
            start = [std(MRS_struct.out.MLalign.ph_aligned.x(ii,:))/2, median(MRS_struct.out.MLalign.ph_aligned.x(ii,:))];
            MRS_struct.out.MLalign.ph_aligned.p(ii,:) = ...
                mle(MRS_struct.out.MLalign.ph_aligned.x(ii,:), 'pdf', Cauchy, 'start', start, 'lower', lb, 'upper', ub, 'options', mleopts);
            MRS_struct.out.MLalign.ph_aligned.fx(ii,:) = ...
                linspace(1.1*min(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 1.1*max(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 1e3);
            MRS_struct.out.MLalign.ph_aligned.pdf(ii,:) = ...
                Cauchy(MRS_struct.out.MLalign.ph_aligned.fx(ii,:), MRS_struct.out.MLalign.ph_aligned.p(ii,1), MRS_struct.out.MLalign.ph_aligned.p(ii,2));
            
            clf(H1);
            subplot(1,2,1);
            bins = linspace(min(MRS_struct.out.MLalign.f_aligned.x(ii,:)), max(MRS_struct.out.MLalign.f_aligned.x(ii,:)), 20);
            binWidth = abs(bins(1)-bins(2));
            h = bar(bins, histcounts(MRS_struct.out.MLalign.f_aligned.x(ii,:), length(bins))/(length(MRS_struct.out.MLalign.f_aligned.x(ii,:))*binWidth), 'histc');
            h.FaceColor = [0.8 0.8 0.8];
            hold on;
            plot(MRS_struct.out.MLalign.f_aligned.fx(ii,:), MRS_struct.out.MLalign.f_aligned.pdf(ii,:), 'Color', [1 0 0], 'LineWidth', 1.2);
            hold off;
            xlabel('\Deltaf (Hz)', 'FontSize', 15);
            ylabel('P(x)', 'FontSize', 15);
            set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
            
            subplot(1,2,2);
            bins = linspace(min(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), max(MRS_struct.out.MLalign.ph_aligned.x(ii,:)), 20);
            binWidth = abs(bins(1)-bins(2));
            h = bar(bins, histcounts(MRS_struct.out.MLalign.ph_aligned.x(ii,:), length(bins))/(length(MRS_struct.out.MLalign.ph_aligned.x(ii,:))*binWidth), 'histc');
            h.FaceColor = [0.8 0.8 0.8];
            hold on;
            plot(MRS_struct.out.MLalign.ph_aligned.fx(ii,:), MRS_struct.out.MLalign.ph_aligned.pdf(ii,:), 'Color', [1 0 0], 'LineWidth', 1.2);
            hold off;
            xlabel('\Delta\phi (deg)', 'FontSize', 15);
            ylabel('P(x)', 'FontSize', 15);
            set(gca, 'FontSize', 12, 'TickDir', 'out', 'Box', 'off');
            
            drawnow;
            %pause(1);
        end
        
        % Line-broadening, zero-filling and FFT
        FullData = DataToAlign .* repmat((exp(-time*MRS_struct.p.LB*pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(FullData, MRS_struct.p.ZeroFillTo(ii), 1),1);
        
        if ~MRS_struct.p.phantom
            
            % Apply a global frequency and zero-order phase shift by fitting a ChoCr model in the frequency domain
            ind = all(MRS_struct.fids.ON_OFF' == 0,2);
            OFF = real(mean(AllFramesFTrealign(:,ind),2));
            
            freqLim  = MRS_struct.spec.freq <= 3.02+0.1 & MRS_struct.spec.freq >= 3.02-0.1;
            [~,i]    = max(abs(OFF(freqLim)));
            freq2    = MRS_struct.spec.freq(freqLim);
            maxFreq  = freq2(i);
            freqLim  = MRS_struct.spec.freq <= maxFreq+0.58 & MRS_struct.spec.freq >= maxFreq-0.42;
            OFF      = OFF(freqLim);
            baseline = (OFF(1) + OFF(end))/2;
            width    = 0.05;
            area     = (max(OFF) - min(OFF)) * width * 4;
            
            x0 = [area width maxFreq 0 baseline 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
            ModelParam = FitChoCr(MRS_struct.spec.freq(freqLim), real(OFF), x0, MRS_struct.p.LarmorFreq(ii));
            
            f   = ModelParam(3) - (3.02 * MRS_struct.p.LarmorFreq(ii));
            phi = ModelParam(4);
            AllFramesFTrealign = ifft(ifftshift(AllFramesFTrealign,1),[],1);
            t = 0:(1/MRS_struct.p.sw(ii)):(size(AllFramesFTrealign,1) - 1) * (1/MRS_struct.p.sw(ii));
            AllFramesFTrealign = AllFramesFTrealign * exp(1i * pi/180 * phi) .* exp(1i * f * 2 * pi .* repmat(t', [1 size(AllFramesFTrealign,2)]));
            AllFramesFTrealign = fftshift(fft(AllFramesFTrealign,[],1),1);
            
            MRS_struct.out.SpecReg.freq{ii}  = MRS_struct.out.SpecReg.freq{ii} + f;
            MRS_struct.out.SpecReg.phase{ii} = MRS_struct.out.SpecReg.phase{ii} + phi;
            
        end
        
        % Reject transients that are greater than 3 st. devs. of MSEs
        % (only applies if not using weighted averaging)
        MRS_struct.out.reject{ii} = zMSE > 3;
        
    end
    
    SpecRegLoop = SpecRegLoop - 1;
    
end

if exist('H1','var')
    close(H1);
end
fprintf('\n');

end



