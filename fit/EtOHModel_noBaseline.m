function [F, J] = EtOHModel_noBaseline(x, freq)
% Function for EtOH model with no baseline
%
% If requested, J returns the analytic Jacobian dF/dx (numel(freq) x 6)

F = x(1) ./ (1 + ((freq - x(2)) / (x(3)/2)).^2) + ...
    x(4) ./ (1 + ((freq - x(5)) / (x(6)/2)).^2);

if nargout > 1
    f = freq(:);
    J = zeros(numel(f), numel(x));
    % Lorentzians: L = A/(1+u^2), u = (f-c)/s, s = FWHM/2
    for ii = [1 4]
        A = x(ii); c = x(ii+1); s = x(ii+2)/2;
        u = (f - c) / s;
        D = 1 + u.^2;
        J(:,ii)   = 1 ./ D;
        J(:,ii+1) = 2 * A * u ./ (s * D.^2);
        J(:,ii+2) = A * u.^2 ./ (s * D.^2); % dL/ds * ds/dFWHM (= 1/2)
    end
end
