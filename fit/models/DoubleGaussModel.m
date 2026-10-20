function F = DoubleGaussModel(x, freq)
% Function for double-Gaussian model

% Two Gaussians
%  x(1) = gaussian amplitude 1
%  x(2) = width 1 ( 1/(2*sigma^2) )
%  x(3) = centre freq of peak 1
%  x(4) = gaussian amplitude 2
%  x(5) = width 2 ( 1/(2*sigma^2) )
%  x(6) = centre freq of peak 2
%  x(7) = amplitude of linear baseline
%  x(8) = constant amplitude offset

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
    x(4) * exp(x(5) * (freq - x(6)).^2) + ...
    x(7) * (freq - x(3)) + x(8);
