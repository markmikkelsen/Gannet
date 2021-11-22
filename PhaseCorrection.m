function fids = PhaseCorrection(fids, MRS_struct)

ii        = MRS_struct.ii;
freqRange = MRS_struct.p.sw(ii)/MRS_struct.p.LarmorFreq(ii);
freq      = (MRS_struct.p.npoints(ii) + 1 - (1:MRS_struct.p.npoints(ii))) / MRS_struct.p.npoints(ii) * freqRange + 4.68 - freqRange/2;
waterLim  = freq <= 4.68 + 0.25 & freq >= 4.68 - 0.25;

spec = fftshift(fft(fids,[],1),1);
if MRS_struct.p.HERMES
    ind = all(MRS_struct.fids.ON_OFF' == 0,2);
else
    switch MRS_struct.p.target{1}
        case {'GABAGlx','GABA','Glx','Lac','EtOH'}
            ind = 1:size(MRS_struct.fids.data,2);
        case 'GSH'
            ind = MRS_struct.fids.ON_OFF == 0;
    end
end
q = sum(abs(spec(waterLim,ind))) * abs(freq(1) - freq(2));
q = (q - median(q)) / median(q) * 100;
q = sum(abs(q) > 40) / length(q);
q_threshold = 0.1;
if q > q_threshold
    water_flag = 1;
else
    water_flag = 0;
end

if MRS_struct.p.HERMES
    n = 4;
else
    n = 2;
end

D      = zeros(size(fids,2)/n);
w      = cell(1,n);
data   = complex(zeros(size(fids,1),n));
time   = (0:(MRS_struct.p.npoints(ii)-1))'/MRS_struct.p.sw(ii);
tMax   = find(time <= 0.1,1,'last');
MSEfun = @(a,b) sum((a - b).^2) / length(a);

if size(MRS_struct.fids.data,2) <= 4
    data = fids;
else
    for jj = 1:n
        if MRS_struct.p.HERMES
            ind = jj:n:size(fids,2);
        else
            ind = find(MRS_struct.fids.ON_OFF == abs(jj-2));
        end
        for kk = 1:size(fids,2)/n
            for ll = 1:size(fids,2)/n
                D(kk,ll) = MSEfun(real(fids(1:tMax,ind(kk))), real(fids(1:tMax,ind(ll))));
            end
        end
        D(~D) = NaN;
        d     = mean(D,'omitnan');
        w{jj} = 1./d.^2;
        w{jj} = w{jj}/sum(w{jj});
        w{jj} = repmat(w{jj}, [size(fids,1) 1]);
        if water_flag
            dataLim = ceil(length(ind)/3);
            data(:,jj) = sum(w{jj}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
        else
            data(:,jj) = sum(w{jj} .* fids(:,ind),2);
        end
    end
end

SUM       = mean(real(fftshift(fft(data,[],1),1)),2);
freqLim   = freq <= 3.02+0.1 & freq >= 3.02-0.1;
[~,i]     = max(abs(SUM(freqLim)));
freq2     = freq(freqLim);
maxFreq   = freq2(i);
freqLim   = freq <= maxFreq+0.58 & freq >= maxFreq-0.42;
SUM_ChoCr = SUM(freqLim);

% Make sure signals have positive phase
[~,i] = max(abs(SUM_ChoCr));
if sign(SUM_ChoCr(i)) < 0
    fids = -fids;
    data = -data;
    SUM  = mean(real(fftshift(fft(data,[],1),1)),2);
    SUM_ChoCr = SUM(freqLim);
end

% Apply zero-order phase correction
Baseline   = (SUM_ChoCr(1) + SUM_ChoCr(end))/2;
Width      = 0.05;
Area       = (max(SUM_ChoCr) - min(SUM_ChoCr)) * Width * 4;
x0         = [Area Width maxFreq 0 Baseline 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
ModelParam = FitChoCr(freq(freqLim), SUM_ChoCr, x0, MRS_struct.p.LarmorFreq(ii));
fids       = fids * exp(1i*pi/180*ModelParam(4));

if size(MRS_struct.fids.data,2) <= 4
    data = fids;
else
    for jj = 1:n
        if MRS_struct.p.HERMES
            ind = jj:4:size(fids,2);
        else
            ind = find(MRS_struct.fids.ON_OFF == abs(jj-2));
        end
        if water_flag
            dataLim = ceil(length(ind)/3);
            data(:,jj) = sum(w{jj}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
        else
            data(:,jj) = sum(w{jj} .* fids(:,ind),2);
        end
    end
end

% Make sure again signals have positive phase
SUM = mean(real(fftshift(fft(data,[],1),1)),2);
SUM_ChoCr = SUM(freqLim);
SUM_ChoCr = SUM_ChoCr - (SUM_ChoCr(1) + SUM_ChoCr(end))/2;
[~,i] = max(abs(SUM_ChoCr));
if sign(SUM_ChoCr(i)) < 0
    fids = -fids;
    data = -data;
    SUM  = mean(real(fftshift(fft(data,[],1),1)),2);
    SUM_ChoCr = SUM(freqLim);
end

% Run phase correction again to make sure the phase is correct
Baseline   = (SUM_ChoCr(1) + SUM_ChoCr(end))/2;
Area       = (max(SUM_ChoCr) - min(SUM_ChoCr)) * Width * 4;
x0         = [Area Width maxFreq 0 Baseline 0 1] .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
ModelParam = FitChoCr(freq(freqLim), SUM_ChoCr, x0, MRS_struct.p.LarmorFreq(ii));
fids       = fids * exp(1i*pi/180*ModelParam(4));

end



