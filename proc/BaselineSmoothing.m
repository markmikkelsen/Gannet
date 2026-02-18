function z = BaselineSmoothing(ii, freq, spec, base_mask, lambda, tol)
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

st = rng; % save current rng
rng('default'); % make sure the legacy rng is not being used
rng(ii); % reproduce the same pseudorandom numbers each time that are unique for each spectrum

if nargin < 6 || isempty(tol)
    tol = 1e-3;
end
if nargin < 5 || isempty(lambda)
    lambda = 1e9;
end

y = spec(:);
base_mask = base_mask(:);
max_iter = 200;
iter = 1;
k = 1;

% Replace non-baseline signal with pseudorandom noise
noise_ind = freq > 9 & freq < 10;
noise_mu  = mean(y(noise_ind));
noise_sd  = std(detrend(y(noise_ind),2)); % detrend noise signal using second-degree polynomial
y(base_mask ~= 1) = noise_mu + noise_sd * randn(sum(base_mask ~= 1),1);
% y(base_mask ~= 1) = 0;

rng(st); % restore previous rng

N = length(y);
D = diff(speye(N), 2); % second-order difference matrix (penalty)
H = lambda * (D' * D);
w = ones(N,1);
% w(base_mask ~= 1) = 0; % set weights to zero for non-baseline signal
alpha = ones(N,1);

show_plots = 0;

if show_plots
    figure(33);
    clf(33);
end

while true    
    A = alpha .* H; % adjust penalty using data-driven coefficient alpha
    W = spdiags(w, 0, N, N);
    
    C = decomposition(W + A, 'banded');
    % We use banded decomposition instead of Cholesky factorization because
    % the penalty matrix won't necessarily be symmetric positive definite
    % and 'banded' is a more efficient solver for banded matrices

    z  = C \ (w .* y);

    d = y - z;

    if show_plots
        ax1 = subplot(4,1,1);
        cla(ax1);
        plot(freq, spec ./ max(spec), 'k', freq, z, 'r');
        set(gca, 'XDir', 'reverse', 'XLim', [-2 7], 'YLim', [-0.25 1.25]);
        drawnow;
    end

    % Get d-, mu(d-), and sigma(d-)
    dn = d(d < 0);
    % m  = mean(dn);
    s  = std(dn);

    % Generalized logistic function of d for weighting
    % wt = 1 ./ (1 + exp(2 * (d - (-m + 2*s)) / s));
    wt = 1 ./ (1 + exp(k * (d - s) / s));
    [wt_sort, ind] = sort(wt);

    if show_plots
        ax2 = subplot(4,1,2);
        cla(ax2);
        hold on;
        % plot(d);
        % plot(w);
        % plot((d - (-m + 2*s)) / s, wt);
        % p = scatter((d - (-m + 2*s)) / s, wt, 'MarkerEdgeColor', '#000', 'MarkerFaceColor', '#1b9e77');
        % p = scatter((d - s) / s, wt, 'MarkerEdgeColor', '#000', 'MarkerFaceColor', '#1b9e77');
        % p.SizeData = 25;
        % plot(wt);
        plot((d(ind) - s) / s, wt_sort, 'LineWidth', 1);
        % plot(dn);
        % plot(d);
        % plot(z);
        plot([1 1], [-0.1 1.1], 'r', 'LineStyle', '--');
        hold off;
        xlim([-4 6]);
        ylim([-0.1 1.1]);
        axis square;
        drawnow;
        xlabel('d');
        ylabel('wt');
        % pause(0.01);

        ax3 = subplot(4,1,3);
        cla(ax3);
        title('wt');
        plot(freq, wt);
        set(gca, 'XDir', 'reverse', 'XLim', [min(freq) max(freq)]);
        drawnow;
        xlabel('ppm');
        ylabel('wt');
    end

    % Check stopping conditions
    if norm(w - wt) / norm(w) < tol || iter == max_iter
        break
    end
    % Update the weights and alpha for the next iteration
    w = wt;
    % w(base_mask ~= 1) = 0;
    alpha = abs(d) / max(abs(d)); % update alpha based on the current residuals
    
    if show_plots
        ax4 = subplot(4,1,4);
        cla(ax4);
        plot(freq, alpha);
        set(gca, 'XDir', 'reverse', 'XLim', [min(freq) max(freq)]);
        drawnow;
        xlabel('ppm');
        ylabel('\alpha');
    end

    iter = iter + 1;
end

end
