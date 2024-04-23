function MRS_struct = SignalAveraging(MRS_struct, AllFramesFT, AllFramesFTrealign, ii, kk, vox)

% Initialize some variables/functions
MSEfun     = @(a,b) sum((a - b).^2) / length(a);
experiment = {'A','B','C','D'};
method     = 'MSE'; % Options: 'MSE', 'MSE2', 'WACFM'

if MRS_struct.p.HERMES
    n = 4;
else
    n = 2;
end

if strcmp(MRS_struct.p.alignment, 'none')
    fprintf('\n');
end

if MRS_struct.p.weighted_averaging && size(MRS_struct.fids.data,2) >= 4 % weighted averaging

    fprintf('Averaging subspectra using weighted averaging and performing subtraction...');
    MRS_struct.p.weighted_averaging_method = method;
    MRS_struct.out.signal_averaging.w{ii}  = zeros(1,size(MRS_struct.fids.data,2));

    freqRange = MRS_struct.p.sw(ii) / MRS_struct.p.LarmorFreq(ii);
    freq = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
    freqLim = freq <= 3.4 & freq >= 1.8;

    for jj = 1:n

        if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
            ind = MRS_struct.fids.ON_OFF == abs(jj-2);
        else
            ind = jj:n:size(AllFramesFTrealign,2);
        end

        % Undo zerofill
        spec = ifft(ifftshift(AllFramesFTrealign(:,ind),1),[],1);
        spec = fftshift(fft(spec(1:MRS_struct.p.npoints(ii),:),[],1),1);

        switch method
            case 'MSE'
                D = zeros(size(AllFramesFTrealign,2)/n);
                for ll = 1:size(AllFramesFTrealign,2)/n
                    for mm = 1:size(AllFramesFTrealign,2)/n
                        D(ll,mm) = MSEfun(real(spec(freqLim,ll)), real(spec(freqLim,mm)));
                    end
                end
                D(~D) = NaN;
                d = median(D,'omitnan');
                w = d.^-2 / sum(d.^-2);
            case 'MSE2'
                d = MSEfun(real(spec(freqLim,:)), median(real(spec(freqLim,:)),2));
                w = d.^-2 / sum(d.^-2);
            case 'WACFM'
                [~,w] = WACFM(real(spec(freqLim,:)), 'GCD');
%                 close(23);
            otherwise
                error('Weighted averaging method not recognized!');
        end
        MRS_struct.out.signal_averaging.w{ii}(ind) = w;
        w = repmat(w, [size(AllFramesFTrealign,1) 1]);
        MRS_struct.spec.(vox{kk}).subspec.(experiment{jj})(ii,:) = sum(w .* AllFramesFTrealign(:,ind),2);

    end

else % conventional averaging

    fprintf('Averaging subspectra and performing subtraction...');
    MRS_struct.p.weighted_averaging = 0; % in case there are 4 or less averages but weighted averaging was still set

    for jj = 1:n
        if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
            ind = MRS_struct.fids.ON_OFF == abs(jj-2);
        else
            ind = jj:n:size(AllFramesFTrealign,2);
        end
        ind = ismember(1:size(AllFramesFTrealign,2), ind);
        if iscell(MRS_struct.out.reject)
            MRS_struct.spec.(vox{kk}).subspec.(experiment{jj})(ii,:) = mean(AllFramesFTrealign(:,ind & MRS_struct.out.reject{ii} == 0),2);
        else
            MRS_struct.spec.(vox{kk}).subspec.(experiment{jj})(ii,:) = mean(AllFramesFTrealign(:,ind & MRS_struct.out.reject(ii) == 0),2);
        end
    end

end

for jj = 1:length(MRS_struct.p.target)

    if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
        ON_ind  = 1;
        OFF_ind = 2;
    else
        ON_ind  = find(MRS_struct.fids.ON_OFF(jj,1:n) == 1);
        OFF_ind = find(MRS_struct.fids.ON_OFF(jj,1:n) == 0);
    end

    if MRS_struct.p.HERMES
        % ON
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:) = ...
            (MRS_struct.spec.(vox{kk}).subspec.(experiment{ON_ind(1)})(ii,:) + ...
            MRS_struct.spec.(vox{kk}).subspec.(experiment{ON_ind(2)})(ii,:)) / 2;
        % OFF
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:) = ...
            (MRS_struct.spec.(vox{kk}).subspec.(experiment{OFF_ind(1)})(ii,:) + ...
            MRS_struct.spec.(vox{kk}).subspec.(experiment{OFF_ind(2)})(ii,:)) / 2;
        % OFF_OFF
        OFF_OFF_ind = all(MRS_struct.fids.ON_OFF(:,1:n)' == 0,2);
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off_off(ii,:) = ...
            MRS_struct.spec.(vox{kk}).subspec.(experiment{OFF_OFF_ind})(ii,:);
    else
        % ON
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:) = ...
            MRS_struct.spec.(vox{kk}).subspec.(experiment{ON_ind})(ii,:);
        % OFF
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:) = ...
            MRS_struct.spec.(vox{kk}).subspec.(experiment{OFF_ind})(ii,:);
    end

    % DIFF
    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:) = ...
        (MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:) - ...
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:)) / 2;

    % DIFF (unaligned)
    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:) = ...
        (mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:) == 1),2) - ...
        mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:) == 0),2)) / 2;

end

end


function [v, w] = WACFM(x, costFun)
% Weighted averaging based on criterion function minimization (WACFM).
% Algorithm from Pander T. A new approach to robust, weighted signal
% averaging. Biocybern Biomed Eng. 2015;35(4):317-327. doi:10.1016/j.bbe.2015.06.002

[~,N]  = size(x);
kStop  = 100;
v      = median(x,2);
w      = zeros(kStop,N);
e      = 1e-6;
m      = 2;
p      = 0.2;
sigma  = 1;
const  = max(abs(v)); % constant to increase robustness to local minima
                      % (﻿Kotowski et al. Biocybern Biomed Eng. 2019. doi:10.1016/j.bbe.2019.09.002)

for k = 1:kStop
    z = x - v;
    w(k,:) = weights(z, costFun, const);

%     figure(23);
%     cla;
%     plot(w(1:k,:)');
%     drawnow;
%     pause(0.25);

    if (k > 1 && norm(w(k,:) - w(k-1,:)) < e) || k == kStop
        w = w(k,:);
        switch costFun % normalize optimal weights so they sum to unity
            case 'square'
                w = w.^m ./ sum(w.^m);
            case 'GCD'
                w = w.^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1 ./ ...
                    sum(w.^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1);
        end
        break
    end

    switch costFun
        case 'square'
            v = sum(w(k,:).^m .* x,2) ./ sum(w(k,:).^m);
        case 'GCD'
            v = sum(w(k,:).^m .* x .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1,2) ./ ...
                sum(w(k,:).^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1);
    end
end

    function w = weights(z, costFun, const)
        switch costFun
            case 'square'
                w = (vecnorm(z) + const).^(2./(1-m)) ./ ...
                    sum((vecnorm(z) + const).^(2./(1-m)));
            case 'GCD'
                w = sum(log(1 + (abs(z) ./ sigma).^p) + const).^(1./(1-m)) ./ ...
                    sum(sum(log(1 + (abs(z) ./ sigma).^p) + const).^(1./(1-m)));
        end
    end

end



