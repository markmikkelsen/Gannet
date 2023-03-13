function F = LorentzGaussModelP(x, freq)
% Function for Lorentz-Gaussian model with phase

% Lorentzian Model multiplied by a Gaussian
% x(1) = Amplitude of (scaled) Lorentzian
% x(2) = 1 / hwhm of Lorentzian (hwhm = half width at half max)
% x(3) = centre freq of Lorentzian
% x(4) = linear baseline slope
% x(5) = constant baseline amplitude
% x(6) =  -1 / 2 * sigma^2  of gaussian
% x(7) = phase (in rad)

% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% F is a normalised Lorentzian - height independent of hwhm
%   = Lorentzian / Peak

F = ((cos(x(7)) * x(1) + sin(x(7)) * x(1) * x(2) * (freq - x(3))) ./ ...
    (x(2)^2 * (freq - x(3)).^2 + 1)) .* ... % Lorentzian
    (exp(x(6) * (freq - x(3)).^2)) + ... % Gaussian
    x(4) * (freq - x(3)) + ... % linear baseline
    x(5); % constant baseline
