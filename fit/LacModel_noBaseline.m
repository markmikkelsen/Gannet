function [F, J] = LacModel_noBaseline(x, freq)
% Function for Lac model with no baseline

% Lac+ (Lorentzians)
%   x(1,4) = amplitudes
%   x(2,5) = widths
%   x(3,6) = freq offsets
% BHB+ (Gaussian)
%   x(7)   = amplitude
%   x(8)   = width
%   x(9)   = freq offset
% MM1.43 (Gaussian)
%   x(10)  = amplitude
%   x(11)  = width
%   x(9)   = freq offset (+0.21 ppm)
%
% If requested, J returns the analytic Jacobian dF/dx (numel(freq) x 11)

% Two Lorentzians + two Gaussians
F = x(1) ./ (1 + ((freq - x(3)) ./ x(2)).^2) + ... % Lac+ (1)
    x(4) ./ (1 + ((freq - x(6)) ./ x(5)).^2) + ... % Lac+ (2)
    x(7) * exp(x(8) * (freq - x(9)).^2) + ... % BHB+
    x(10) * exp(x(11) * (freq - (x(9) + 0.21)).^2); % MM1.43

if nargout > 1
    f = freq(:);
    J = zeros(numel(f), numel(x));
    % Lorentzians: L = A/(1+u^2), u = (f-c)/s
    for ii = [1 4]
        A = x(ii); s = x(ii+1); c = x(ii+2);
        u = (f - c) / s;
        D = 1 + u.^2;
        J(:,ii)   = 1 ./ D;
        J(:,ii+1) = 2 * A * u.^2 ./ (s * D.^2);
        J(:,ii+2) = 2 * A * u ./ (s * D.^2);
    end
    % Gaussians (shared center x(9))
    [~, J(:,7), J(:,8), dc1]   = GaussTermJacobian(x(7),  x(8),  x(9),        f);
    [~, J(:,10), J(:,11), dc2] = GaussTermJacobian(x(10), x(11), x(9) + 0.21, f);
    J(:,9) = dc1 + dc2;
end
