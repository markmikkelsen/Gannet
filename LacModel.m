function F = LacModel(x,freq)

% Lac+
%   x(1) = gaussian amplitude
%   x(2) = 1/(2*sigma^2)
%   x(3) = centre freq of peak
% BHB+
%   x(4) = gaussian amplitude
%   x(5) = 1/(2*sigma^2)
%   x(6) = centre freq of peak
% MM3
%   x(10 = gaussian amplitude
% Baseline
%   x(7) = offset
%   x(8) = slope
%   x(9) = quadratic

F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
    x(1) * exp(x(2) * (freq - (x(3) - 0.055)).^2) + ...
    x(4) * exp(x(5) * (freq - x(6)).^2) + ...    
    x(10) * exp(x(5) * (freq - (x(6) + 0.21)).^2) + ...
    x(7) + x(8) * (freq - x(3)) + x(9) * (freq - x(3)).^2;

% F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
%     x(1) * exp(x(2) * (freq - (x(3) - 0.055)).^2) + ...
%     x(4) * exp(x(5) * (freq - (x(6))).^2) + ...    
%     x(7) + x(8) * (freq - x(6)) + x(9) * (freq - x(6)).^2;



