function [baseline, signal] = BaselineRecognition(y, freq)

% Power spectrum of first-derivative of signal calculated by CWT
y = abs(cwt_ricker(y, 75)).^2;

noiseLim = freq <= 9 & freq >= 8;
sigma    = std(y(noiseLim));

w = 1:5;
k = 3;
[baseline, signal] = deal(zeros(size(y)));

while 1
    if w(end) > length(y)
        break
    end
    h = max(y(w)) - min(y(w));
    if h < k*sigma
        baseline(w) = 1;
    else
        signal(w) = 1;
    end
    w = w + 1;
end

end


function [W, scales, t] = cwt_ricker(x, scales, dt, nSigma)
% CWT_RICKER  Continuous wavelet transform (CWT) of 1-D data using Ricker wavelet.
%
%   [W, scales, t] = cwt_ricker(x, scales, dt, nSigma)
%
% Inputs
%   x      : 1-D signal (vector)
%   scales : vector of positive scales (e.g., logspace(log10(2),log10(128),40))
%   dt     : sample spacing (default = 1)
%   nSigma : support half-width in std devs (default = 5). Larger = more accurate, slower.
%
% Outputs
%   W      : CWT coefficients, size [numel(scales) x numel(x)]
%   scales : the scales used (returned for convenience)
%   t      : time axis (same length as x), in units of dt
%
% Notes
% - Ricker mother wavelet:
%       psi(u) = (2/(sqrt(3)*pi^(1/4))) * (1 - u^2) * exp(-u^2/2)
%   (unit-energy normalization in continuous time).
% - CWT definition used:
%       W(a,b) = (1/sqrt(a)) * ∫ x(t) * psi((t-b)/a) dt
%   Implemented via discrete convolution with dt scaling.
%
% Example
%   x = chirp(0:999,0,999,0.2) + 0.2*randn(1,1000);
%   scales = logspace(log10(2), log10(128), 40);
%   [W,sc,t] = cwt_ricker(x, scales, 1);
%   imagesc(t, sc, abs(W)); axis xy; xlabel('t'); ylabel('scale'); colorbar;
%
% MM (251218): Created using ChatGPT 5.2

if nargin < 3 || isempty(dt)
    dt = 1;
end

if nargin < 4 || isempty(nSigma)
    nSigma = 5;
end

x = x(:).';
N = numel(x);
t = (0:N-1) * dt;

if nargin < 2 || isempty(scales)
    scales = logspace(log10(1), log10(max(2, floor(N/8))), 40);
end
scales = scales(:);
if any(scales <= 0)
    error('All scales must be positive.');
end

% Mother wavelet normalization constant (unit energy)
C = 2/(sqrt(3)*pi^(1/4));

W = zeros(numel(scales), N);

for j = 1:numel(scales)
    a = scales(j);

    % Build a time grid for the wavelet support: u = (tau)/a
    % Choose tau range = ±(nSigma * a) (since exp(-u^2/2) decays fast)
    tauMax = nSigma * a;       % in same units as t
    tau = (-tauMax:dt:tauMax); % lag axis
    u = tau / a;

    psi = C * (1 - u.^2) .* exp(-0.5*u.^2); % mother psi(u)

    % Scale the wavelet: psi_a(tau) = (1/sqrt(a)) * psi(tau/a)
    psi_a = (1/sqrt(a)) * psi;

    % Discrete approx to integral adds dt
    psi_a = psi_a * dt;

    % Correlation via convolution with flipped wavelet
    W(j,:) = conv(x, fliplr(psi_a), 'same');
end

end
