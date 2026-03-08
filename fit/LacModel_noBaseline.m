function F = LacModel_noBaseline(x, freq)
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
%   x(9)   = freq offset

% Two Lorentzians + two Gaussians
F = x(1) ./ (1 + ((freq - x(3)) ./ x(2)).^2) + ... % Lac+ (1)
    x(4) ./ (1 + ((freq - x(6)) ./ x(5)).^2) + ... % Lac+ (2)
    x(7) * exp(x(8) * (freq - x(9)).^2) + ... % BHB+
    x(10) * exp(x(11) * (freq - (x(9) + 0.21)).^2); % MM1.43
