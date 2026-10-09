function [z, lambda_opt, rmse_e] = BaselineSmoothing(freq, spec, lambda, tol, ext, rule)
% Estimate a smoothed baseline with extended range based on penalized least
% squares (erPLS), an extension of adaptive smoothness parameter penalized
% least squares (asPLS) and asymmetrically reweighted penalized least
% squares (arPLS)
%
% Baek et al. Baseline correction using asymmetrically reweighted penalized
%   least squares smoothing. Analyst. 2015;140(1):250-257.
%   doi:10.1039/C4AN01061B
% Zhang et al. Baseline correction for infrared spectra using adaptive
%   smoothness parameter penalized least squares method. Spectrosc Lett.
%   2020;53(3):222-233. doi:10.1080/00387010.2020.1730908
% Zhang et al. An automatic baseline correction method based on the
%   penalized least squares method. Sensors. 2020;20(7):2015.
%   doi:10.3390/s20072015
%
% Usage:
%   z = BaselineSmoothing(freq, spec, lambda, tol)
%       asPLS baseline with a fixed smoothing parameter lambda (default:
%       1e9)
%
%   [z, lambda_opt, rmse_e] = BaselineSmoothing(freq, spec, lambda, tol, ext, rule)
%       erPLS baseline. ext is the output of ExtendSpectrum(freq, spec),
%       which carries out steps 1-3 (linear extension of one end of the
%       spectrum, with noise, and addition of positive and negative test
%       peaks). Here, lambda is a vector of candidate values (default: a
%       0.1-decade grid from about 10^-6 to 10^5 times fwhm_pts^4, where fwhm_pts
%       is the test-peak FWHM in points). For each candidate, asPLS is run
%       on the windowed, extended spectrum and the RMSE between the fitted
%       baseline and the noise-free linear extension is computed in the
%       extended range (RMSE_e; step 4). lambda is then selected (step 5)
%       and used to estimate the baseline of the full original spectrum
%       (step 6). rmse_e is returned for every candidate (NaN for
%       candidates skipped by the coarse-to-fine search)
%
%       rule sets how lambda is selected in step 5:
%         'knee' (default) - the smallest lambda for which RMSE_e has
%                            dropped below the noise SD of the spectrum
%                            edge (ext.noise_sd), i.e., the most flexible
%                            baseline that no longer follows the test
%                            peaks. Used for MRS because the flat ends of
%                            MRS spectra make RMSE_e keep falling as lambda
%                            grows, so its minimum always picks the
%                            stiffest candidate
%         'min'            - the lambda with the lowest RMSE_e, as in
%                            Zhang et al.

if nargin < 4 || isempty(tol)
    tol = 1e-4;
end

if nargin < 5 || isempty(ext)
    % asPLS with a fixed smoothing parameter
    if nargin < 3 || isempty(lambda)
        lambda = 1e9;
    end
    z          = asPLS(freq, spec, lambda, tol);
    lambda_opt = lambda;
    rmse_e     = [];
    return
end

% erPLS (steps 4-6)
if nargin < 6 || isempty(rule)
    rule = 'knee';
end
rule = validatestring(rule, {'knee', 'min'});
if nargin < 3 || isempty(lambda)
    % lambda scales roughly with (width in points)^4, so center the grid on
    % the test-peak FWHM
    c = 4 * log10(ext.fwhm_pts);
    lambda = 10.^((floor(c) - 6):0.1:(ceil(c) + 5));
end
lambda = sort(lambda(:));
nl     = length(lambda);
if length(spec) ~= ext.N
    error('ext must be the output of ExtendSpectrum for the same spectrum.');
end

% Step 4: Calculate RMSE_e for each candidate lambda. To save time, a
% coarse subset of the candidates is searched first, then every candidate
% between the coarse neighbors of the coarse selection
rmse_e = NaN(nl,1);
stride = 5;
if nl > 2*stride
    coarse = unique([1:stride:nl, nl]);
else
    coarse = 1:nl;
end
for ll = coarse
    rmse_e(ll) = ExtendedRangeRMSE(ext, lambda(ll), tol);
end
[ind, found] = SelectLambda(rmse_e, rule, ext.noise_sd);
if length(coarse) < nl
    pos  = find(coarse == ind);
    fine = coarse(max(pos-1, 1)):coarse(min(pos+1, length(coarse)));
    fine = fine(isnan(rmse_e(fine)));
    for ll = fine
        rmse_e(ll) = ExtendedRangeRMSE(ext, lambda(ll), tol);
    end
    [ind, found] = SelectLambda(rmse_e, rule, ext.noise_sd);
end

% Step 5: Select lambda
lambda_opt = lambda(ind);
if ~found
    warning('BaselineSmoothing:noKnee', ...
        'RMSE_e never dropped below the noise SD; using the lambda with the lowest RMSE_e (%.3g).', lambda_opt);
elseif nl > 1 && (ind == 1 || ind == nl)
    warning('BaselineSmoothing:lambdaAtBound', ...
        'Selected lambda (%.3g) is at the edge of the search range; consider widening it.', lambda_opt);
end

% Step 6: Estimate the baseline of the original spectrum with the optimal lambda
z = asPLS(freq, spec, lambda_opt, tol);

end


function [ind, found] = SelectLambda(rmse_e, rule, noise_sd)
% Index of the selected lambda among the candidates evaluated so far (NaN
% entries are ignored). found is false if the 'knee' rule had no crossing
[~, ind_min] = min(rmse_e);
found = true;
if strcmp(rule, 'min')
    ind = ind_min;
    return
end
% Last candidate at or below the minimum whose RMSE_e exceeds the noise SD;
% the next evaluated candidate is the selection. Restricting the search to
% lambdas up to the minimum ignores the erratic, ill-conditioned fits at the
% very stiff end of the grid
evaluated = find(~isnan(rmse_e));
evaluated = evaluated(evaluated <= ind_min);
above     = evaluated(rmse_e(evaluated) > noise_sd);
if isempty(above)
    ind = evaluated(1);
elseif above(end) == ind_min
    ind   = ind_min; % never dropped below the noise SD
    found = false;
else
    ind = evaluated(find(evaluated > above(end), 1));
end
end


function r = ExtendedRangeRMSE(ext, lambda, tol)
% RMSE between the asPLS baseline of the extended spectrum and the
% noise-free linear extension (without the test peaks) in the extended
% range (Eq. 5)
z = asPLS(ext.freq, ext.spec, lambda, tol);
r = sqrt(mean((ext.line - z(ext.ind)).^2));
end


function z = asPLS(freq, spec, lambda, tol)
% asPLS baseline estimate with a fixed smoothing parameter lambda

y         = spec(:);
max_iter  = 400;
iter      = 1;
k         = 0.5;
w_min     = 1e-6;                    % weight floor; keeps W + A non-singular
s_min     = 1e3 * eps(max(abs(y)));  % scale-aware floor on sigma(d-)

N = length(y);
D = diff(speye(N), 2); % second-order difference matrix (penalty)
H = lambda * (D' * D);
w = ones(N,1);
alpha = ones(N,1);

show_plots = 0;

if show_plots
    close all;
    figure(33);
end

while true    
    A = alpha .* H; % adjust penalty using data-driven coefficient alpha
    W = spdiags(w, 0, N, N);    
    C = decomposition(W + A, 'banded');
    % We use banded decomposition instead of Cholesky factorization because
    % the penalty matrix won't necessarily be symmetric positive definite
    % and 'banded' is a more efficient solver for banded matrices
    z = C \ (w .* y);
    
    d = y - z;

    if show_plots
        ax2 = subplot(4,1,2);
        cla(ax2);
        hold on;
        plot(freq, spec ./ max(spec), 'k');
        plot(freq, y ./ max(spec), 'b');
        plot(freq, z ./ max(spec), 'r');
        hold off;
        set(gca, 'XDir', 'reverse', 'XLim', [-2 7], 'YLim', [-0.5 1.25]);
        xlabel('ppm');
        drawnow;
    end

    % Iteration cap reached: z above was refitted with the weights and alpha
    % from the final update, so return it without computing new ones
    if iter > max_iter
        break
    end

    % Get d-, mu(d-), and sigma(d-)
    dn = d(d < 0);
    if numel(dn) < 2
        % Baseline lies on or below the data everywhere, so sigma(d-) is
        % undefined and there is nothing left to push down; keep current z
        break
    end
    % m = mean(dn);
    s = std(dn);
    if ~isfinite(s) || s < s_min
        s = s_min;
    end

    % Generalized logistic function of d for weighting
    % wt = 1 ./ (1 + exp(2 * (d - (-m + 2*s)) / s));
    % Written as 1/(1+exp(x)) = (1 - tanh(x/2))/2, which cannot overflow when
    % d >> s. The floor stops a weight of exactly 0 coinciding with alpha ~ 0
    % and zeroing a row of W + A (singular system)
    wt = 0.5 * (1 - tanh(k * (d - s) / (2*s)));
    wt = max(wt, w_min);
    [wt_sort, ind] = sort(wt);

    if show_plots
        ax1 = subplot(4,1,1);
        cla(ax1);
        hold on;
        plot((d(ind) - s) / s, wt_sort, 'LineWidth', 1);
        plot([1 1], [-0.1 1.1], 'r', 'LineStyle', '--');
        hold off;
        xlim([-4 6]);
        ylim([-0.1 1.1]);
        axis square;
        xlabel('d');
        ylabel('wt');
        drawnow;

        ax3 = subplot(4,1,3);
        cla(ax3);
        title('wt');
        plot(freq, wt);
        set(gca, 'XDir', 'reverse', 'XLim', [-2 7]);
        xlabel('ppm');
        ylabel('wt');
        drawnow;
    end

    % Check convergence. If it isn't reached within max_iter updates, the loop
    % runs one more solve with the last weights (see the iter > max_iter check
    % above) so the returned z is consistent with them, unlike standard arPLS
    if norm(w - wt) / norm(w) < tol
        break
    end
    % Update the weights and alpha for the next iteration
    w = wt;
    d_max = max(abs(d));
    if d_max > 0
        alpha = abs(d) / d_max; % update alpha based on the current residuals
    else
        alpha = ones(N,1);      % exact fit: fall back to uniform smoothness
    end

    if show_plots
        ax4 = subplot(4,1,4);
        cla(ax4);
        plot(freq, alpha);
        set(gca, 'XDir', 'reverse', 'XLim', [-2 7]);
        xlabel('ppm');
        ylabel('\alpha');
        drawnow;
    end

    iter = iter + 1;
end

end
