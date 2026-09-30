function [F, J] = GABAGlxModel_noBaseline(x, freq)
% Function for GABA+Glx model with no baseline

%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = center freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = center freq of peak 2
%  x(7) = gaussian amplitude 3
%  x(8) = width 3 ( 1/(2*sigma^2) )
%  x(9) = center freq of peak 3
%
% If requested, J returns the analytic Jacobian dF/dx (numel(freq) x 9)

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
    x(4) * exp(x(5) * (freq - x(6)).^2) + ...
    x(7) * exp(x(8) * (freq - x(9)).^2);

if nargout > 1
    f = freq(:);
    J = zeros(numel(f), numel(x));
    for k = 1:3
        ii = 3*(k-1) + (1:3);
        [~, J(:,ii(1)), J(:,ii(2)), J(:,ii(3))] = GaussTermJacobian(x(ii(1)), x(ii(2)), x(ii(3)), f);
    end
end
