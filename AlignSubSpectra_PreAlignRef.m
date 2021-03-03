function MRS_struct = AlignSubSpectra_PreAlignRef(MRS_struct, water_flag)

ii        = MRS_struct.ii;
fids      = MRS_struct.fids.data_align;
fids_ref  = MRS_struct.fids.data;
freqRange = MRS_struct.p.sw(ii) / MRS_struct.p.LarmorFreq(ii);
freq      = (size(fids,1) + 1 - (1:size(fids,1))) / size(fids,1) * freqRange + 4.68 - freqRange/2;
t         = 0:(1/MRS_struct.p.sw(ii)):(size(fids,1)-1)*(1/MRS_struct.p.sw(ii));
tMax      = find(t <= 0.1,1,'last');
MSEfun    = @(a,b) sum((a - b).^2) / length(a);

% Pre-allocate memory
if MRS_struct.p.HERMES
    n = 4;
else
    n = 2;
end
D        = zeros(size(fids,2)/n);
D_ref    = zeros(size(fids,2)/n);
w        = cell(1,n);
w_ref    = cell(1,n);
data     = complex(zeros(size(fids,1),n));
data_ref = complex(zeros(size(fids,1),n));

% Optimization options
lsqnonlinopts = optimoptions(@lsqnonlin);
lsqnonlinopts = optimoptions(lsqnonlinopts,'Algorithm','levenberg-marquardt','Display','off');

% Use weighted averaging to average subspectra
if size(MRS_struct.fids.data,2) <= 4
    data     = fids;
    data_ref = fids_ref;
else
    for jj = 1:n
        
        if MRS_struct.p.HERMES
            ind = jj:n:size(fids,2);
        else
            ind = find(MRS_struct.fids.ON_OFF == abs(jj-2));
        end
        
        for kk = 1:size(fids,2)/n
            for ll = 1:size(fids,2)/n
                D(kk,ll)     = feval(MSEfun, real(fids(1:tMax,ind(kk))), real(fids(1:tMax,ind(ll))));
                D_ref(kk,ll) = feval(MSEfun, real(fids_ref(1:tMax,ind(kk))), real(fids_ref(1:tMax,ind(ll))));
            end
        end
        
        D(~D) = NaN;
        d     = nanmean(D);
        w{jj} = 1./d.^2;
        w{jj} = w{jj} / sum(w{jj});
        w{jj} = repmat(w{jj}, [size(fids,1) 1]);
        
        D_ref(~D_ref) = NaN;
        d_ref         = nanmean(D_ref);
        w_ref{jj}     = 1./d_ref.^2;
        w_ref{jj}     = w_ref{jj} / sum(w_ref{jj});
        w_ref{jj}     = repmat(w_ref{jj}, [size(fids,1) 1]);
        
        if water_flag
            dataLim        = ceil(length(ind)/3);
            data(:,jj)     = sum(w{jj}(:,1:dataLim) .* fids(:,ind(1:dataLim)),2);
            data_ref(:,jj) = sum(w_ref{jj}(:,1:dataLim) .* fids_ref(:,ind(1:dataLim)),2);
        else
            data(:,jj)     = sum(w{jj} .* fids(:,ind),2);
            data_ref(:,jj) = sum(w_ref{jj} .* fids_ref(:,ind),2);
        end
        
    end
end

flatdata     = FlattenData(data);
ref_flatdata = FlattenData(data_ref);

% Set HERMES subexperiment indices (A, B, C, D)
if MRS_struct.p.HERMES
    if ~MRS_struct.p.HERCULES
        if length(MRS_struct.p.target) == 2 && all(strcmp(MRS_struct.p.target,{'GABAGlx','GSH'}))
            switch MRS_struct.p.vendor
                case 'GE'
                    subSpecInd = [3 2 1 4];
                case {'Philips','Philips_data','Philips_raw'}
                    subSpecInd = [1 2 3 4];
                case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                    subSpecInd = [3 1 4 2];
            end
        elseif length(MRS_struct.p.target) == 3 && all(strcmp(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            switch MRS_struct.p.vendor
                case {'Philips','Philips_data','Philips_raw'}
                    % throw an error for now
                case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                    subSpecInd = [3 1 4 2];
            end
        end
    else
        switch MRS_struct.p.vendor
            case {'Philips','Philips_data','Philips_raw'}
                subSpecInd = [1 4 3 2];
            case {'Siemens_twix','Siemens_rda','Siemens_dicom'}
                subSpecInd = [3 2 1 4];
        end
    end
end

a = max(flatdata(:));
freqLim = freq <= 3.4 & freq >= 1.8;

% Align subspectra to their pre-aligned equivalent
if MRS_struct.p.HERMES
    
    % A
    in         = flatdata(:,:,subSpecInd([1 4]))./a;
    in(:,:,2)  = ref_flatdata(:,:,subSpecInd(1))./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(1,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
    % B
    in         = flatdata(:,:,subSpecInd([2 4]))./a;
    in(:,:,2)  = ref_flatdata(:,:,subSpecInd(2))./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(2,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
    % C
    in         = flatdata(:,:,subSpecInd([3 4]))./a;
    in(:,:,2)  = ref_flatdata(:,:,subSpecInd(3))./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(3,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
    % D
    in         = flatdata(:,:,subSpecInd([4 1]))./a;
    in(:,:,2)  = ref_flatdata(:,:,subSpecInd(4))./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(4,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
else
    
    % ON
    in(:,:,1)  = flatdata(:,:,1)./a;
    in(:,:,2)  = ref_flatdata(:,:,1)./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(1,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
    % OFF
    in(:,:,1)  = flatdata(:,:,2)./a;
    in(:,:,2)  = ref_flatdata(:,:,2)./a;
    fun        = @(x) objFunc(in, freq, freqLim, t, x);
    param(2,:) = lsqnonlin(fun, [0 0], [], [], lsqnonlinopts);
    
end

% Apply frequency and phase corrections to subspectra
if MRS_struct.p.HERMES
    
    ind1 = subSpecInd(1):4:size(fids,2);
    ind2 = subSpecInd(2):4:size(fids,2);
    ind3 = subSpecInd(3):4:size(fids,2);
    ind4 = subSpecInd(4):4:size(fids,2);
    
    fids(:,ind1) = fids(:,ind1) ...
                   .* exp(1i * param(1,1) * 2 * pi .* repmat(t', [1 length(ind1)])) ...
                   .* repmat(exp(1i * pi/180 * param(1,2)), [length(t) length(ind1)]);
    fids(:,ind2) = fids(:,ind2) ...
                   .* exp(1i * param(2,1) * 2 * pi .* repmat(t', [1 length(ind2)])) ...
                   .* repmat(exp(1i * pi/180 * param(2,2)), [length(t) length(ind2)]);
    fids(:,ind3) = fids(:,ind3) ...
                   .* exp(1i * param(3,1) * 2 * pi .* repmat(t', [1 length(ind3)])) ...
                   .* repmat(exp(1i * pi/180 * param(3,2)), [length(t) length(ind3)]);
    fids(:,ind4) = fids(:,ind4) ...
                   .* exp(1i * param(4,1) * 2 * pi .* repmat(t', [1 length(ind4)])) ...
                   .* repmat(exp(1i * pi/180 * param(4,2)), [length(t) length(ind4)]);
    
    MRS_struct.out.SpecReg.freq{ii}(ind1)  = MRS_struct.out.SpecReg.freq{ii}(ind1) + param(1,1);
    MRS_struct.out.SpecReg.phase{ii}(ind1) = MRS_struct.out.SpecReg.phase{ii}(ind1) + param(1,2);
    MRS_struct.out.SpecReg.freq{ii}(ind2)  = MRS_struct.out.SpecReg.freq{ii}(ind2) + param(2,1);
    MRS_struct.out.SpecReg.phase{ii}(ind2) = MRS_struct.out.SpecReg.phase{ii}(ind2) + param(2,2);
    MRS_struct.out.SpecReg.freq{ii}(ind3)  = MRS_struct.out.SpecReg.freq{ii}(ind3) + param(3,1);
    MRS_struct.out.SpecReg.phase{ii}(ind3) = MRS_struct.out.SpecReg.phase{ii}(ind3) + param(3,2);
    MRS_struct.out.SpecReg.freq{ii}(ind4)  = MRS_struct.out.SpecReg.freq{ii}(ind4) + param(4,1);
    MRS_struct.out.SpecReg.phase{ii}(ind4) = MRS_struct.out.SpecReg.phase{ii}(ind4) + param(4,2);
    
else
    
    ind1 = find(MRS_struct.fids.ON_OFF == 1);
    ind2 = find(MRS_struct.fids.ON_OFF == 0);    
    
    fids(:,ind1) = fids(:,ind1) ...
                   .* exp(1i * param(1,1) * 2 * pi .* repmat(t', [1 length(ind1)])) ...
                   .* repmat(exp(1i * pi/180 * param(1,2)), [length(t) length(ind1)]);
    fids(:,ind2) = fids(:,ind2) ...
                   .* exp(1i * param(2,1) * 2 * pi .* repmat(t', [1 length(ind2)])) ...
                   .* repmat(exp(1i * pi/180 * param(2,2)), [length(t) length(ind2)]);
    
    MRS_struct.out.SpecReg.freq{ii}(ind1)  = MRS_struct.out.SpecReg.freq{ii}(ind1) + param(1,1);
    MRS_struct.out.SpecReg.phase{ii}(ind1) = MRS_struct.out.SpecReg.phase{ii}(ind1) + param(1,2);
    MRS_struct.out.SpecReg.freq{ii}(ind2)  = MRS_struct.out.SpecReg.freq{ii}(ind2) + param(2,1);
    MRS_struct.out.SpecReg.phase{ii}(ind2) = MRS_struct.out.SpecReg.phase{ii}(ind2) + param(2,2);
    
end

MRS_struct.fids.data_align = fids;

if ishandle(44)
    close(44);
end

end


function [flatdata, data] = FlattenData(data)

flatdata(:,1,:) = real(data);
flatdata(:,2,:) = imag(data);
data = real(fftshift(fft(data,[],1),1));

end


function out = objFunc(in, freq, freqLim, t, x) %#ok<INUSL>

f   = x(1);
phi = x(2);

y = complex(in(:,1,1), in(:,2,1));
y = y .* exp(1i * pi * (t' * f * 2 + phi/180));

a = real(fftshift(fft(y)));
b = real(fftshift(fft(complex(in(:,1,2), in(:,2,2)))));

DIFF = a - b;
out  = DIFF(freqLim);

% figure(44);
% cla;
% hold on;
% plot(freq, a, 'k');
% plot(freq, b, 'r');
% plot(freq, DIFF - 5, 'k');
% plot(freq(freqLim), out - 5, 'r');
% hold off;
% set(gca,'XDir','reverse','XLim',[1 5]);
% drawnow;
% pause(0.05);

end



