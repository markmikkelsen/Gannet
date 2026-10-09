function ext = ExtendSpectrum(freq, spec, opts)
% Linearly extend one end of a spectrum and add a Gaussian peak to the
% extended range. This implements steps 1-3 of the extended range penalized
% least squares (erPLS) method. The output is passed to BaselineSmoothing,
% which carries out steps 4-6 (selection of the optimal smoothing parameter
% lambda of asPLS)
%
%   Step 1: Fit a first-order polynomial to the last Omega points at the
%           chosen end of the spectrum (Omega = N/20 by default)
%   Step 2: Extrapolate the line over W new points (W = N/5 by default) to
%           give the extended signal y_e
%   Step 3: Add a Gaussian peak y_g with height H (max(spec) by default)
%           and width W/2, centered in the extended range, to give
%           y_eg = y_e + y_g. The width is taken as the full +/-3 sigma
%           span (sigma = W/12), so the outer W/4 on each side of the
%           peak is (practically) the bare line, which anchors the
%           baseline at the far end of the extended range. The new
%           spectrum is y_new = [y, y_eg]
%
% Inputs:
%   freq - frequency axis (ppm); must be uniformly spaced
%   spec - real-valued spectrum (same length as freq)
%   Optional name-value arguments:
%     Side        - end of the spectrum to extend: 'upfield' (default) or
%                   'downfield'
%     FitFraction - length of the linear fitting range (Omega) as a fraction
%                   of N (default: 1/20)
%     ExtFraction - length of the extended range (W) as a fraction of N
%                   (default: 1/5)
%     Height      - height of the added Gaussian peak (H) (default:
%                   max(spec))
%
% Output (structure):
%   ext.freq  - extended frequency axis (column vector)
%   ext.spec  - extended spectrum, y_new (column vector)
%   ext.line  - linear extension without the Gaussian (y_e); the reference
%               for the extended-range RMSE (RMSE_e)
%   ext.gauss - added Gaussian peak (y_g)
%   ext.ind   - indices of the extended range in ext.freq/ext.spec
%   ext.orig  - indices of the original spectrum in ext.freq/ext.spec
%
% Zhang et al. An automatic baseline correction method based on the
%   penalized least squares method. Sensors. 2020;20(7):2015.
%   doi:10.3390/s20072015

arguments
    freq {mustBeVector, mustBeReal, mustBeFinite}
    spec {mustBeVector, mustBeReal, mustBeFinite}
    opts.Side {mustBeMember(opts.Side, {'upfield', 'downfield'})} = 'upfield'
    opts.FitFraction (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.FitFraction, 1)} = 1/20
    opts.ExtFraction (1,1) double {mustBePositive} = 1/5
    opts.Height (1,1) double {mustBeReal} = NaN
end

freq = double(freq(:));
y    = double(spec(:));
N    = length(y);

if length(freq) ~= N
    error('freq and spec must have the same length.');
end

n_fit = max(2, round(opts.FitFraction * N)); % length of Omega
W     = max(3, round(opts.ExtFraction * N)); % length of extended range
df    = (freq(end) - freq(1)) / (N - 1);     % signed frequency step

% Gannet frequency axes run from downfield to upfield, but don't assume it
upfield_at_end = freq(end) < freq(1);
extend_at_end  = strcmp(opts.Side, 'upfield') == upfield_at_end;

if extend_at_end
    fit_ind = (N - n_fit + 1):N;
    freq_e  = freq(end) + df * (1:W).';
else
    fit_ind = 1:n_fit;
    freq_e  = freq(1) - df * (W:-1:1).';
end

% Step 1: Linear fit over Omega (centered and scaled for conditioning)
[p, ~, mu] = polyfit(freq(fit_ind), y(fit_ind), 1);

% Step 2: Linear expansion
y_e = polyval(p, freq_e, [], mu);

% Step 3: Signal addition
if isnan(opts.Height)
    H = max(y);
    if H <= 0
        % Spectrum is entirely non-positive; use its largest magnitude so
        % the added peak still points upward, as asPLS assumes
        H = max(abs(y));
    end
else
    H = opts.Height;
end
sigma = (W / 2) / 6; % +/-3 sigma spans W/2
x     = (1:W).' - (W + 1)/2;
y_g   = H * exp(-x.^2 / (2 * sigma^2));
y_eg  = y_e + y_g;

if extend_at_end
    ext.freq = [freq; freq_e];
    ext.spec = [y; y_eg];
    ext.ind  = (N + 1:N + W).';
    ext.orig = (1:N).';
else
    ext.freq = [freq_e; freq];
    ext.spec = [y_eg; y];
    ext.ind  = (1:W).';
    ext.orig = (W + 1:W + N).';
end
ext.line  = y_e;
ext.gauss = y_g;

end
