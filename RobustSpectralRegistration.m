function [AllFramesFTrealign, MRS_struct] = RobustSpectralRegistration(MRS_struct)
% Align using robust spectral registration (Mikkelsen et al. NMR Biomed.
% 2020;33(10):e4368. doi:10.1002/nbm.4368)

ii = MRS_struct.ii;

% Looping parameters
if MRS_struct.p.HERMES
    SpecRegLoop = 3;
    SubspecToAlign = repmat([3 2 1 0], [1 size(MRS_struct.fids.data,2)/4]);
else
    SpecRegLoop = 1;
    SubspecToAlign = MRS_struct.fids.ON_OFF;
end

% Pre-allocate memory
if MRS_struct.p.HERMES
    params = zeros(size(MRS_struct.fids.data,2)/4,2);
    MSE    = zeros(1,size(MRS_struct.fids.data,2)/4);
    w      = zeros(1,size(MRS_struct.fids.data,2)/4);
else
    params = zeros(size(MRS_struct.fids.data,2)/2,2);
    MSE    = zeros(1,size(MRS_struct.fids.data,2)/2);
    w      = zeros(1,size(MRS_struct.fids.data,2)/2);
end
MRS_struct.out.SpecReg.freq{ii}  = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.phase{ii} = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.out.SpecReg.MSE{ii}   = zeros(1,size(MRS_struct.fids.data,2));
zMSE = zeros(1,size(MRS_struct.fids.data,2));
MRS_struct.fids.data_align = complex(zeros(size(MRS_struct.fids.data)));
DataToAlign = complex(zeros(size(MRS_struct.fids.data)));
MSEfun = @(a,b) sum((a - b).^2) / length(a);

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');

% Automatic unstable lipid/residual water removal
freqRange = MRS_struct.p.sw(ii) / MRS_struct.p.LarmorFreq(ii);
freq      = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
waterLim  = freq <= 4.68+0.25 & freq >= 4.68-0.25;
lipidLim  = freq <= 1.85 & freq >= 0;
noiseLim  = freq <= 9 & freq >= 8;

spec = fftshift(fft(MRS_struct.fids.data,[],1),1);

S     = mean(real(spec),2);
noise = S(noiseLim);
fit   = polyval(polyfit(freq(noiseLim), noise.', 2), freq(noiseLim)); % detrend noise
noise = noise - fit.';
r     = std(S(lipidLim)) / std(noise);

if MRS_struct.p.HERMES
    ind  = all(MRS_struct.fids.ON_OFF' == 0,2);
    ind2 = ind;
else
    switch MRS_struct.p.target{1}
        case {'GABAGlx','GABA','Glx','Lac','EtOH'}
            ind  = 1:size(MRS_struct.fids.data,2);
            ind2 = MRS_struct.fids.ON_OFF == 0;
        case 'GSH'
            ind  = MRS_struct.fids.ON_OFF == 0;
            ind2 = ind;
    end
end

q = sum(abs(spec(waterLim,ind))) * abs(freq(1) - freq(2));
q = (q - median(q)) / median(q) * 100;
q = sum(abs(q) > 40) / length(q);

% Force-run water filter if very strong water suppression was used
S        = mean(abs(spec(:,ind2)),2);
maxNAA   = max(S(freq <= 4.25 & freq >= 1.8));
maxWater = max(S(freq <= 4.68+0.22 & freq >= 4.68-0.22));
if maxWater / maxNAA < 1.5
    q = 1;
end

r_threshold = 40;
q_threshold = 0.1;
lipid_flag  = 0;
water_flag  = 0;

if r > r_threshold || q > q_threshold
    if r > r_threshold
        lipid_flag = 1;
    end
    if q > q_threshold
        water_flag = 1;
    end
    
    % Turn off warnings about the legacy random number generator if they are currently on
    w1 = warning('query','MATLAB:RandStream:ActivatingLegacyGenerators');
    w2 = warning('query','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
    if strcmp(w1.state,'on')
        warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
    end
    if strcmp(w2.state,'on')
        warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
    end
    
    reverseStr = '';
    for jj = 1:size(MRS_struct.fids.data,2)
        if lipid_flag && ~water_flag
            msg = sprintf('\nUnstable lipid contamination detected. Temporarily removing lipids from transient: %d', jj);
        elseif ~lipid_flag && water_flag
            msg = sprintf('\nUnstable residual water detected. Temporarily removing residual water from transient: %d', jj);
        elseif lipid_flag && water_flag
            msg = sprintf('\nUnstable lipid contamination and residual water detected.\nTemporarily removing lipids and residual water from transient: %d', jj);
        end
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        DataToAlign(:,jj) = SignalFilter(spec(:,jj), lipid_flag, water_flag, jj, MRS_struct);
    end
    
    % Turn warnings back on if they were previously on
    if strcmp(w1.state,'on')
        warning('on','MATLAB:RandStream:ActivatingLegacyGenerators');
    end
    if strcmp(w2.state,'on')
        warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
    end
else
    DataToAlign = MRS_struct.fids.data;
end

time = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);

% Spectral registration
count = 1;
reverseStr = '';
while SpecRegLoop > -1
    
    % Use first n points of time-domain data, where n is the last point where abs(diff(mean(SNR))) > 0.5
    signal = abs(DataToAlign(:,SubspecToAlign == SpecRegLoop));
    noise  = 2*std(signal(ceil(0.75*size(signal,1)):end,:));
    SNR    = signal ./ repmat(noise, [size(DataToAlign,1) 1]);
    SNR    = abs(diff(mean(SNR,2)));
    SNR    = SNR(time <= 0.2); % use no more than 200 ms of data
    tMax   = find(SNR > 0.5,1,'last');
    if isempty(tMax) || tMax < find(time <= 0.05,1,'last') % use at least 50 ms of data
                                                           % (shortened this from 100 ms because seems
                                                           % like this helps when there are spurious echoes
                                                           % or strong water suppression was used)
        tMax = find(time <= 0.05,1,'last');
    end
    
    % Flatten complex data for use in spectral registration
    clear flatdata
    flatdata(:,1,:) = real(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    flatdata(:,2,:) = imag(DataToAlign(1:tMax,SubspecToAlign == SpecRegLoop));
    
    % Determine optimal alignment order by calculating a similarity metric (mean squared error)
    if size(MRS_struct.fids.data,2) <= 4
        alignOrder = 1;
    else
        D = zeros(size(flatdata,3));
        ind = find(SubspecToAlign == SpecRegLoop);
        for jj = 1:size(flatdata,3)
            for kk = 1:size(flatdata,3)
                D(jj,kk) = feval(MSEfun, real(DataToAlign(1:tMax,ind(jj))), real(DataToAlign(1:tMax,ind(kk))));
            end
        end
        D(~D) = NaN;
        d = median(D,'omitnan');
        [~,alignOrder] = sort(d);
    end
    
    % Set initial reference transient based on similarity index
    target = squeeze(flatdata(:,:,alignOrder(1)));
    target = target(:);
    
    % Scalar to normalize transients (reduces optimization time)
    a = max(abs(target));
    
    % Pre-allocate memory
    m = zeros(length(target), size(flatdata,3));
    
    % Starting values for optimization
    f0   = MRS_struct.spec.F0freq2{ii}(ind) * MRS_struct.p.LarmorFreq(ii);
    f0   = f0(alignOrder);
    f0   = f0 - f0(1);
    phi0 = zeros(size(f0));
    x0   = [f0(:) phi0(:)];
    
    % Determine frequency and phase offsets by spectral registration
    t = 0:(1/MRS_struct.p.sw(ii)):(length(target)/2-1)*(1/MRS_struct.p.sw(ii));
    iter = 1;
    for jj = alignOrder
        msg = sprintf('\nRunning robust spectral registration on transient: %d', count);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
        count = count + 1;
        
        transient    = squeeze(flatdata(:,:,jj));
        fun          = @(x) SpecReg(transient(:)/a, target/a, t, x);
        params(jj,:) = lsqnonlin(fun, x0(iter,:), [], [], lsqnonlinopts);
        
        f       = params(jj,1);
        phi     = params(jj,2);
        m_c     = complex(flatdata(:,1,jj), flatdata(:,2,jj));
        m_c     = m_c .* exp(1i*pi*(t'*f*2+phi/180));
        m(:,jj) = [real(m_c); imag(m_c)];
        resid   = target - m(:,jj);
        MSE(jj) = sum(resid.^2) / (length(resid) - 2);
        
        % Update reference
        w(jj)  = 0.5*corr(target, m(:,jj)).^2;
        target = (1 - w(jj))*target + w(jj)*m(:,jj);
        
        iter = iter + 1;
    end
    
    ind = find(SubspecToAlign == SpecRegLoop);
    MRS_struct.out.SpecReg.freq{ii}(ind)  = params(:,1);
    MRS_struct.out.SpecReg.phase{ii}(ind) = params(:,2);
    MRS_struct.out.SpecReg.MSE{ii}(ind)   = MSE;
    zMSE(ind(MSE > 0)) = zscore(MSE(MSE > 0)); % standardized MSEs (exclude the first reference)
    
    % Apply frequency and phase corrections to raw data
    MRS_struct.fids.data_align(:,ind) = MRS_struct.fids.data(:,ind) ...
                                        .* exp(repmat(1i * params(:,1)' * 2 * pi, [length(time) 1]) .* repmat(time, [1 size(flatdata,3)])) ...
                                        .* repmat(exp(1i * pi/180 * params(:,2)'), [length(time) 1]);
    
    if SpecRegLoop == 0
        
        % Align subspectra
        if MRS_struct.p.use_prealign_ref
            MRS_struct = AlignSubSpectra_PreAlignRef(MRS_struct, water_flag);
        else
            MRS_struct = AlignSubSpectra(MRS_struct, water_flag, r);
        end
        
        % Line-broadening, zero-filling and FFT
        AllFramesFTrealign = MRS_struct.fids.data_align .* repmat((exp(-time * MRS_struct.p.LB * pi)), [1 size(MRS_struct.fids.data,2)]);
        AllFramesFTrealign = fftshift(fft(AllFramesFTrealign, MRS_struct.p.ZeroFillTo(ii), 1),1);
        
        % Global frequency shift
        % Unclear if this is now redundant as it's already done in
        % AlignSubSpectra above; leaving alone for now (MM: 210330)
        if ~MRS_struct.p.phantom
            CrFreqRange        = MRS_struct.spec.freq <= 3.02+0.15 & MRS_struct.spec.freq >= 3.02-0.15;
            [~,FrameMaxPos]    = max(abs(mean(real(AllFramesFTrealign(CrFreqRange,:)),2)));
            freq               = MRS_struct.spec.freq(CrFreqRange);
            CrFreqShift        = freq(FrameMaxPos);
            CrFreqShift        = CrFreqShift - 3.02;
            CrFreqShift_pts    = round(CrFreqShift / abs(MRS_struct.spec.freq(1) - MRS_struct.spec.freq(2)));
            AllFramesFTrealign = circshift(AllFramesFTrealign, CrFreqShift_pts, 1);
            MRS_struct.out.SpecReg.freq{ii} = MRS_struct.out.SpecReg.freq{ii} + CrFreqShift * MRS_struct.p.LarmorFreq(ii);
        end
        
        % Reject transients that are greater than 3 st. devs. of MSEs
        % (only applies if not using weighted averaging)
        MRS_struct.out.reject{ii} = zMSE > 3;
        
    end
    
    SpecRegLoop = SpecRegLoop - 1;
    
end

fprintf('\n');

end



