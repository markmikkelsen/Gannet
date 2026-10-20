function F = LorentzGaussModel(x, freq)
% Function for Lorentz-Gaussian model

% CJE 24Nov10 - removed phase term from fit - this is now dealt with
% by the phasing of the water ref scans in MRSLoadPfiles
% Lorentzian Model multiplied by a Gaussian.
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = (x(1) ./ (x(2)^2 * (freq - x(3)).^2 + 1)) ... % Lorentzian
    .* (exp(x(6) * (freq - x(3)).^2)) ... % Gaussian
    + x(4) * (freq - x(3)) ... % linear baseline
    + x(5); % constant baseline
