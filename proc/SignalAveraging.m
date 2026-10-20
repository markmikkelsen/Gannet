function MRS_struct = SignalAveraging(MRS_struct, AllFramesFT, AllFramesFTrealign, ii, kk, vox)

% Initialize some variables/functions
MSEfun = @(a,b) sum((a - b).^2) / length(a);
method = 'MSE'; % Options: 'MSE', 'MSE2', 'WACFM'

if MRS_struct.p.HERMES
    n_subspec = 4;
    experiment = {'A','B','C','D'};
else
    n_subspec = 2;
    experiment = {'ON','OFF'};
    if ~MRS_struct.fids.ON_OFF(1) % OFF is first
        experiment = flip(experiment);
    end
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

    for jj = 1:n_subspec

        if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
            ind = MRS_struct.fids.ON_OFF == abs(jj-2);
        else
            ind = jj:n_subspec:size(AllFramesFTrealign,2);
        end

        % Undo zerofill
        spec = ifft(ifftshift(AllFramesFTrealign(:,ind),1),[],1);
        spec = fftshift(fft(spec(1:MRS_struct.p.npoints(ii),:),[],1),1);

        switch method
            case 'MSE'
                D = zeros(size(AllFramesFTrealign,2)/n_subspec);
                for ll = 1:size(AllFramesFTrealign,2)/n_subspec
                    for mm = 1:size(AllFramesFTrealign,2)/n_subspec
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
                % close(23);
            otherwise
                error('Weighted averaging method not recognized!');
        end
        MRS_struct.out.signal_averaging.w{ii}(ind) = w;
        w = repmat(w, [size(AllFramesFTrealign,1) 1]);
        MRS_struct.spec.(vox{kk}).subspec.(experiment{jj})(ii,:) = sum(w .* AllFramesFTrealign(:,ind),2);

    end

else % arithmetic averaging

    fprintf('Averaging subspectra using arithmetic averaging and performing subtraction...');
    MRS_struct.p.weighted_averaging = 0; % in case there are 4 or less averages but weighted averaging was still set

    for jj = 1:n_subspec
        if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
            ind = MRS_struct.fids.ON_OFF == abs(jj-2);
        else
            ind = jj:n_subspec:size(AllFramesFTrealign,2);
        end
        ind = ismember(1:size(AllFramesFTrealign,2), ind);
        MRS_struct.spec.(vox{kk}).subspec.(experiment{jj})(ii,:) = mean(AllFramesFTrealign(:,ind & MRS_struct.out.reject{ii} == 0),2);
    end

end

for jj = 1:length(MRS_struct.p.target)

    if strcmp(MRS_struct.p.vendor, 'Philips') && strcmp(MRS_struct.p.seqorig, 'Philips')
        ON_ind  = 1;
        OFF_ind = 2;
    else
        ON_ind  = find(MRS_struct.fids.ON_OFF(jj,1:n_subspec) == 1);
        OFF_ind = find(MRS_struct.fids.ON_OFF(jj,1:n_subspec) == 0);
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
        OFF_OFF_ind = all(MRS_struct.fids.ON_OFF(:,1:n_subspec)' == 0,2);
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
    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_unfilt(ii,:) = ...
        (MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:) - ...
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:)) / 2;

    % DIFF (unaligned)
    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_unfilt_noalign(ii,:) = ...
        (mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:) == 1),2) - ...
        mean(AllFramesFT(:,MRS_struct.fids.ON_OFF(jj,:) == 0),2)) / 2;

    % SUM
    MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).sum_unfilt(ii,:) = ...
        (MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).on(ii,:) + ...
        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).off(ii,:)) / 2;

    % Plot averaged subspectra
    if MRS_struct.p.plot_subspec
        PlotAvgSubpec(MRS_struct, ii, kk, vox, n_subspec, experiment, jj)
    end

end

end
