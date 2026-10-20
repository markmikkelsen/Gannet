function [F, J] = EightGaussModel_noBaseline(x, freq)
% Function for eight-Gaussian model with no baseline

% x(1)     = Gaussian 1 amplitude
% x(2)     = Gaussian 1 width (1/(2*sigma^2))
% x(3)     = Gaussian 1 center freq
% x(4-5)   = Gaussian 2
% x(6-8)   = Gaussian 3
% x(9-10)  = Gaussian 4
% x(11-13) = Gaussian 5
% x(14-15) = Gaussian 6
% x(16-18) = Gaussian 7
% x(19-21) = Gaussian 8
%
% If requested, J returns the analytic Jacobian dF/dx (numel(freq) x 21)

F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
    x(4) * exp(x(5) * (freq - (x(8) + 0.1)).^2) + ...
    x(6) * exp(x(7) * (freq - x(8)).^2) + ...
    x(9) * exp(x(10) * (freq - (x(13) + 0.05)).^2) + ...
    x(11) * exp(x(12) * (freq - x(13)).^2) + ...
    x(14) * exp(x(15) * (freq - (x(13) - 0.05)).^2) + ...
    x(16) * exp(x(17) * (freq - x(18)).^2) + ...
    x(19) * exp(x(20) * (freq - x(21)).^2);

if nargout > 1
    f = freq(:);
    J = zeros(numel(f), numel(x));
    % Each Gaussian: [amplitude, width, center] parameter indices and
    % center offset (Gaussians 2, 4 and 6 share centers with 3 and 5)
    idx = [ 1  2  3;
            4  5  8;
            6  7  8;
            9 10 13;
           11 12 13;
           14 15 13;
           16 17 18;
           19 20 21];
    offset = [0 0.1 0 0.05 0 -0.05 0 0];
    for k = 1:size(idx,1)
        [~, dA, dw, dc] = GaussTermJacobian(x(idx(k,1)), x(idx(k,2)), x(idx(k,3)) + offset(k), f);
        J(:,idx(k,1)) = dA;
        J(:,idx(k,2)) = dw;
        J(:,idx(k,3)) = J(:,idx(k,3)) + dc; % accumulate for shared centers
    end
end
