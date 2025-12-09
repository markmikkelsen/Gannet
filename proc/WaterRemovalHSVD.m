function FIDsvd = WaterRemovalHSVD(FID, sw, nHSVD, HSVDLow, HSVDHigh, plottingon, usepoints)
%***********************************************************
% WaterRemovalHSVD.m
%
% Residual water resonance removal by SVD.
%
% 1. SVD of spectral region [ZxLow ... ZxHigh]
% 2. Least-squares fitting to model the total spectrum
% 3. Extraction of resonances in the region [WxLow ... WxHigh]
% 4. Subtraction of water from original FID
%
% By Robin A. de Graaf
% MRRC, Yale University
% Original : March, 1999
% Modified : September, 2008
%***********************************************************
if plottingon
    disp('starting SVD water removal.....')
end

DataHandling = 3;

zff = 1;

FIDorig = FID;
clear FID;

nporig = length(FIDorig(:));
% npsvd = nporig;
npsvd = usepoints; % use different number of points for speed

FID = FIDorig(1:npsvd);
FIDsvd = FID(1:npsvd);

npnonzero = length(FIDsvd);

tacq = npnonzero/(1000*sw);                % Acquisition time (in s)
dt = tacq/npnonzero;                       % Dwell-time (in s)
time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
time = reshape(time,npnonzero,1);

Lmax = round(0.4*npnonzero);               % Dimension 1 for LxM SVD matrix
Mmax = npnonzero+1-Lmax;                   % Dimension 2 for LxM SVD matrix

if plottingon
    disp(' ');
    disp('Water removal in progress ... ');
end

%*****************
% Allocate memory
%*****************
H = zeros(Lmax,Mmax);

%**********************************************
% Create Hankel matrix from original FID data
%**********************************************
for L = 1:Lmax
    M = 1:Mmax;
    H(L,M) = FIDsvd(L+M-1);
end

%**********************************************
% Perform SVD on Hankel matrix
%**********************************************
if plottingon
    tic;
    disp('Step 1 : SVD of Hankel matrix in progress ...');
end

[U,~,~] = svd(H);

if plottingon
    tt = toc;
    dd = ['... done in ' num2str(tt,3) ' s.'];
    disp(dd);
end

Uup = zeros(Lmax-1,nHSVD);
Udown = zeros(Lmax-1,nHSVD);
% Udownc = zeros(nHSVD,Lmax-1);

if plottingon
    disp('Step 2 : Calculation of lineshape parameters in progress ...');
    tic;
end

%**********************************************
% Calculate truncated SVD matrix
%**********************************************
Utr = U(:,1:nHSVD);
% Str = S(:,1:nHSVD);
% Vtr = V(:,1:nHSVD);

for kk1 = 2:Lmax
    for kk2 = 1:nHSVD
        Uup(kk1-1,kk2) = Utr(kk1,kk2);
        Udown(kk1-1,kk2) = Utr(kk1-1,kk2);
    end
end

Z = pinv(Udown)*Uup;
q = log(eig(Z));

%******************************************************
% Determination of frequencies and T2 constants from D
%******************************************************
% Frequency (in Hz)
frq = imag(q)/(2*pi*dt);
% Time constant (in s)
decay = real(q)/dt;

switch DataHandling
    case 1
        time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
        time = reshape(time,npnonzero,1);
        basis = zeros(length(time),nHSVD);
        
        % Calculate basis functions
        for kk1 = 1:nHSVD
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDsvd;
    case 2
        time = 0:dt:(npsvd-1)*dt;              % Time base of FID
        time = reshape(time,npsvd,1);
        basis = zeros(length(time),nHSVD);
        
        % Calculate basis functions
        for kk1 = 1:nHSVD
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDsvd;
    case 3
        time = 0:dt:(nporig-1)*dt;              % Time base of FID
        time = reshape(time,nporig,1);
        basis = zeros(length(time),nHSVD);
        
        % Calculate basis functions
        for kk1 = 1:nHSVD
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDorig;
end

amp = abs(ampcomplex);
phs = atan2(imag(ampcomplex),real(ampcomplex));

if plottingon
    tt = toc;
    dd1 = [' ... done in ' num2str(tt,3) ' s.'];
    disp(dd1);
    
    disp('Step 3 : Spectral reconstruction in progress ...');
    
    tic;
end

%***********************************************************************
% Reconstruct signal in the frequency range [HSVDLow, HSVDHigh]
%***********************************************************************
waterpos = find(((frq > -1000*HSVDHigh) & (frq < -1000*HSVDLow)));

nwater = length(waterpos);

switch DataHandling
    case 1
        FIDw = zeros(npnonzero,1);
        %FIDw0 = zeros(npnonzero,1);
    case 2
        FIDw = zeros(npsvd,1);
        %FIDw0 = zeros(npsvd,1);
    case 3
        FIDw = zeros(nporig,1);
        %FIDw0 = zeros(nporig,1);
end

% meanfrq = mean(frq(waterpos(1:nwater)));

for kk1 = 1:nwater
    FIDcomponent = amp(waterpos(kk1)).*exp(2*pi*1i*frq(waterpos(kk1))*time).*exp(time*decay(waterpos(kk1))).*exp(1i*phs(waterpos(kk1)));
    %FIDcomponent0 = amp(waterpos(kk1)).*exp(2*pi*1i*(frq(waterpos(kk1))-meanfrq)*time).*exp(time*decay(waterpos(kk1))).*exp(1i*phs(waterpos(kk1)));
    FIDw = FIDw + FIDcomponent;
    %FIDw0 = FIDw0 + FIDcomponent0;
    clear FIDcomponent FIDcomponent0;
end

switch DataHandling
    case 1
        if npsvd ~= npnonzero
            FIDw(npsvd+1:npsvd) = 0;
            %FIDw0(npsvd+1:npsvd) = 0;
            time = 0:dt:(npsvd-1)*dt;              % Time base of FID
            tacq = (npsvd-1)*dt;
        end
    case 2
        if npsvd ~= npnonzero
            FIDw(npsvd+1:nporig) = 0;
            %FIDw0(npsvd+1:nporig) = 0;
            time = 0:dt:(nporig-1)*dt;              % Time base of FID
            tacq = (nporig-1)*dt;
        else
            FIDw(npsvd+1:nporig) = 0.0;
            %FIDw0(npsvd+1:nporig) = 0.0;
            time = 0:dt:(nporig-1)*dt;              % Time base of FID
            tacq = (nporig-1)*dt;
        end
        FID = FIDorig;
    case 3
        FID = FIDorig;
end

% Fourier transformation
spec = fftshift(fft(FIDorig, zff*nporig));
spec = reshape(spec,length(spec),1);
specw = fftshift(fft(FIDw, zff*nporig));
specw = reshape(specw,length(specw),1);

% Phase correction
specA = real(spec);
specAw = real(specw);

if plottingon
    tt = toc;
    dd = ['... done, using ' num2str(nwater) ' water components, in ' num2str(tt,3) ' s.'];
    disp(dd);
    disp('Water removal completed.');
end

%*********************************************************************
% Display original and fitted FID and spectra and the differences
%*********************************************************************
if plottingon
    hh = figure;
    set(hh,'position',[200 50 800 600])
    
    freq = 0.5*sw:(-sw/(zff*nporig-1)):-0.5*sw;
    
    subplot(2,2,1), plot(time,real(FID),'b',time,real(FIDw),'r');
    axis([0 tacq 1.1*min(real(FID)) 1.1*max(real(FID))])
    title('Original (blue)/fitted (red) FID');
    subplot(2,2,3), plot(time,real(FID-FIDw));
    axis([0 tacq 1.1*min(real(FID-FIDw)) 1.1*max(real(FID-FIDw))])
    title('Difference FID');
    
    subplot(2,2,2), plot(freq,specA,'b',freq,specAw,'r');
    %axis([ZxLow ZxHigh ZyLow ZyHigh])
    title('Original (blue)/fitted (red) spectrum');
    subplot(2,2,4), plot(freq,(specA-specAw));
    %axis([ZxLow ZxHigh ZyLow ZyHigh])
    title('Difference spectrum');
end

FIDsvd = FID - FIDw;

end



