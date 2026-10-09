function ext = ExtendSpectrum(freq, spec, opts)
% Extend one end of an MRS spectrum with a straight line plus noise and add
% a pair of synthetic test peaks (one positive, one negative) to the
% extended range. This implements steps 1-3 of the extended range penalized
% least squares (erPLS) method, adapted for in vivo 1H MRS. The output is
% passed to BaselineSmoothing, which carries out steps 4-6 (selection of the
% smoothing parameter lambda of asPLS)
%
% Adaptations for in vivo 1H MRS (Zhang et al. developed erPLS for IR and
% Raman spectra):
%   - All sizes are set in ppm rather than as fractions of the number of
%     points, so the result does not depend on spectral width, zero-filling
%     or field strength
%   - The extension carries a mirrored copy of the spectrum's own edge
%     noise (the residual of the linear fit), so asPLS sees realistic,
%     correctly correlated noise in the extended range
%   - Two test peaks are added: a positive one with the height of the most
%     positive signal in PeakRange and a negative one with the height of the
%     most negative signal there (if any). Edited difference spectra contain
%     large negative signals (e.g., NAA, Asp), which a baseline must not
%     follow either
%   - Only the last WindowWidth ppm of the spectrum are kept next to the
%     extension. The smoothing parameter only needs local data to be
%     calibrated, and the much shorter signal makes the lambda search fast
%
%   Step 1: Fit a first-order polynomial to the last FitWidth ppm at the
%           chosen end of the spectrum
%   Step 2: Extrapolate the line over ExtWidth ppm and add the mirrored
%           fit residual (noise) to give the extended signal y_e
%   Step 3: Add the Gaussian test peaks y_g (FWHM in ppm) at 1/3 and 2/3 of
%           the extended range, giving y_eg = y_e + y_g. The new spectrum is
%           y_new = [y(window), y_eg]
%
% Inputs:
%   freq - frequency axis (ppm); must be uniformly spaced
%   spec - real-valued spectrum (same length as freq)
%   Optional name-value arguments:
%     Side        - end of the spectrum to extend: 'upfield' (default) or
%                   'downfield'
%     FitWidth    - width of the linear fitting range, in ppm (default: 2)
%     ExtWidth    - width of the extended range, in ppm (default: 2; widened
%                   to 10*FWHM if needed so the test peaks stay separated)
%     WindowWidth - width of the original spectrum kept next to the
%                   extension, in ppm (default: 3; Inf keeps all of it)
%     FWHM        - FWHM of the Gaussian test peaks, in ppm (default: 0.15)
%     PeakRange   - ppm range used to set the test-peak heights
%                   (default: [0.5 4.25])
%     Height      - test-peak heights [H_pos H_neg]; a scalar gives a
%                   positive peak only (default: [max min] of the spectrum
%                   in PeakRange; the negative peak is omitted if min >= 0)
%     Noise       - add the mirrored edge noise to the extension (default:
%                   true)
%
% Output (structure):
%   ext.freq     - frequency axis of the windowed and extended spectrum
%   ext.spec     - windowed and extended spectrum, y_new (column vector)
%   ext.line     - noise-free linear extension without the test peaks; the
%                  reference for the extended-range RMSE (RMSE_e)
%   ext.gauss    - added test peaks (y_g)
%   ext.ind      - indices of the extended range in ext.freq/ext.spec
%   ext.win      - indices of the kept window in the original spectrum
%   ext.N        - length of the original spectrum
%   ext.noise_sd - standard deviation of the edge noise (fit residual)
%   ext.fwhm_pts - FWHM of the test peaks, in points
%
% Zhang et al. An automatic baseline correction method based on the
%   penalized least squares method. Sensors. 2020;20(7):2015.
%   doi:10.3390/s20072015

arguments
    freq {mustBeVector, mustBeReal, mustBeFinite}
    spec {mustBeVector, mustBeReal, mustBeFinite}
    opts.Side {mustBeMember(opts.Side, {'upfield', 'downfield'})} = 'upfield'
    opts.FitWidth (1,1) double {mustBePositive} = 2
    opts.ExtWidth (1,1) double {mustBePositive} = 2
    opts.WindowWidth (1,1) double {mustBePositive} = 3
    opts.FWHM (1,1) double {mustBePositive} = 0.15
    opts.PeakRange (1,2) double {mustBeReal} = [0.5 4.25]
    opts.Height double {mustBeReal} = []
    opts.Noise (1,1) logical = true
end

freq = double(freq(:));
y    = double(spec(:));
N    = length(y);

if length(freq) ~= N
    error('freq and spec must have the same length.');
end

df  = (freq(end) - freq(1)) / (N - 1); % signed frequency step
ppp = 1 / abs(df);                     % points per ppm

n_fit = max(2, round(opts.FitWidth * ppp));
n_fit = min(n_fit, N);
W     = round(max(opts.ExtWidth, 10 * opts.FWHM) * ppp);
n_win = min(N, max(n_fit, round(opts.WindowWidth * ppp)));

% Gannet frequency axes run from downfield to upfield, but don't assume it
upfield_at_end = freq(end) < freq(1);
extend_at_end  = strcmp(opts.Side, 'upfield') == upfield_at_end;

if extend_at_end
    fit_ind = (N - n_fit + 1):N;
    win     = (N - n_win + 1):N;
    freq_e  = freq(end) + df * (1:W).';
else
    fit_ind = 1:n_fit;
    win     = 1:n_win;
    freq_e  = freq(1) - df * (W:-1:1).';
end

pk_lims = sort(opts.PeakRange);
if any(freq(fit_ind) >= pk_lims(1) & freq(fit_ind) <= pk_lims(2))
    warning('ExtendSpectrum:fitOverlapsPeaks', ...
        'The linear fitting range overlaps PeakRange; reduce FitWidth or extend the other side.');
end

% Step 1: Linear fit (centered and scaled for conditioning)
[p, ~, mu] = polyfit(freq(fit_ind), y(fit_ind), 1);
res = y(fit_ind) - polyval(p, freq(fit_ind), [], mu);

% Step 2: Linear expansion plus mirrored edge noise. Everything below is
% built outward from the edge of the spectrum and flipped at the end if the
% extension goes before the first point
y_e = polyval(p, freq_e, [], mu);
if extend_at_end
    res_out = flipud(res); % edge point first
else
    res_out = res;
end
if opts.Noise
    noise_out = repmat([res_out; flipud(res_out)], ceil(W / (2*n_fit)), 1);
    noise_out = noise_out(1:W);
else
    noise_out = zeros(W,1);
end

% Step 3: Signal addition
if isempty(opts.Height)
    in_range = freq >= pk_lims(1) & freq <= pk_lims(2);
    if ~any(in_range)
        error('PeakRange does not overlap the frequency axis.');
    end
    H = [max(y(in_range)), min(y(in_range))];
else
    H = opts.Height(:).';
end
H_pos = max(H(1), 0);
if numel(H) > 1
    H_neg = min(H(2), 0);
else
    H_neg = 0;
end
fwhm_pts = opts.FWHM * ppp;
sigma    = fwhm_pts / (2 * sqrt(2 * log(2)));
x        = (1:W).';
y_g_out  = H_pos * exp(-(x - W/3).^2 / (2 * sigma^2)) + ...
           H_neg * exp(-(x - 2*W/3).^2 / (2 * sigma^2));

if extend_at_end
    noise = noise_out;
    y_g   = y_g_out;
    ext.freq = [freq(win); freq_e];
    ext.spec = [y(win); y_e + noise + y_g];
    ext.ind  = (n_win + 1:n_win + W).';
else
    noise = flipud(noise_out);
    y_g   = flipud(y_g_out);
    ext.freq = [freq_e; freq(win)];
    ext.spec = [y_e + noise + y_g; y(win)];
    ext.ind  = (1:W).';
end
ext.line     = y_e;
ext.gauss    = y_g;
ext.win      = win(:);
ext.N        = N;
ext.noise_sd = std(res);
ext.fwhm_pts = fwhm_pts;

end
