function z = BaselineSmoothing(ii, freq, y, b, target, lambda, ratio)
% Estimate a smoothed baseline with asymmetrically reweighted penalized
% least squares (arPLS)
%
% Baek et al. Baseline correction using asymmetrically reweighted penalized
% least squares smoothing. Analyst. 2015;140(1):250-257.
% doi:10.1039/C4AN01061B

s = rng; % save current rng
rng('default'); % make sure the legacy rng is not being used
rng(ii); % reproduce the same pseudorandom numbers each time that are unique for each spectrum

if nargin < 7
    ratio = 1e-5;
    if nargin < 6
        lambda = 1e3;
    end
end

y = y(:);
b = b(:);
k_stop = 500;
k = 1;
new_lambda = false;

% Replace non-baseline signal with pseudorandom noise
noise_mu  = mean(y(freq > 9 & freq < 10));
noise_sd  = std(detrend(y(freq > 9 & freq < 10),2));
y(b ~= 1) = noise_mu + noise_sd * randn(sum(b ~= 1),1);

rng(s); % restore previous rng

% switch target
% 
%     case 'GABAGlx'
%         freqLim.Glx  = [3.52 4.0];
%         freqLim.GABA = [2.77 3.26];
%         freqLim.NAA  = [1.76 2.54];
%         freqLim.MM09 = [0.61 1.18];
% 
%         y(freq < freqLim.Glx(2) & freq > freqLim.Glx(1))   = 0;
%         y(freq < freqLim.GABA(2) & freq > freqLim.GABA(1)) = 0;
%         y(freq < freqLim.NAA(2) & freq > freqLim.NAA(1))   = 0;
%         y(freq < freqLim.MM09(2) & freq > freqLim.MM09(1)) = 0;
% 
%     case 'GSH'
%         freqLim.mI  = [3.39 4.10];
%         freqLim.GSH = [2.15 3.24];
%         freqLim.Lac = [1.05 1.6];
% 
%         y(freq < freqLim.mI(2) & freq > freqLim.mI(1))   = 0;
%         y(freq < freqLim.GSH(2) & freq > freqLim.GSH(1)) = 0;
%         y(freq < freqLim.Lac(2) & freq > freqLim.Lac(1)) = 0;
% 
%     case 'SUM'
%         freqLim.water = [4.32 5.13];
%         freqLim.mI = [3.54 4.10];
%         freqLim.CrCho = [2.9 3.34];
%         freqLim.NAA = [1.76 2.21];
% 
%         y(freq < freqLim.water(2) & freq > freqLim.water(1)) = noise_mu;
%         y(freq < freqLim.mI(2) & freq > freqLim.mI(1)) = noise_mu;
%         y(freq < freqLim.CrCho(2) & freq > freqLim.CrCho(1)) = noise_mu;
%         y(freq < freqLim.NAA(2) & freq > freqLim.NAA(1)) = noise_mu;
% end

N = length(y);
D = diff(speye(N), 2); % second-order difference matrix (penalty)
H = lambda * (D' * D);
w = ones(N,1);
w(b ~= 1) = 0; % set weights to zero for non-baseline signal

% figure(33);
% clf(33);
% % plot(y,'k','LineWidth',1);
% hold on;
% set(gca,'xlim',[3.5e4 5.5e4]);

while true
    if new_lambda
        H = lambda * (D' * D);
        new_lambda = false;
    end

    % Cholesky decomposition
    W = spdiags(w, 0, N, N);
    try
        C = chol(W + H);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:posdef')
            warning('Matrix is not positive definite; decreasing lambda from %g to %g...', lambda, lambda / 10);
            lambda = lambda / 10; % decrease lambda for the next iteration
            new_lambda = true;
            continue
        end
    end

    z = C \ (C' \ (w .* y));
    d = y - z;

    % Get d-, mu(d-), and sigma(d-)
    dn = d(d < 0);
    m  = mean(dn);
    s  = std(dn);

    % Generalized logistic function of d for weighting
    wt = 1 ./ (1 + exp(2 * (d - (-m + 2*s)) / s));

    % clf(33);
    % hold on;
    % plot(d);
    % plot(w);
    % plot(d, wt);
    % plot(sort(wt));
    % plot(dn);
    % plot(d);
    % plot(z);
    % hold off;
    % drawnow;
    % pause(0.5);

    % Check stopping conditions
    if norm(w - wt) / norm(w) < ratio || k > k_stop
        break
    end
    w = wt;
    k = k + 1;
end
