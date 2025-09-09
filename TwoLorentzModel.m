function F = TwoLorentzModel(x, freq)
% Function for double-Lorentzian model

% CJE LorentzModel
% Lorentzian  = (1/pi) * (hwhm) / (deltaf^2 + hwhm^2) (Wolfram)
% Peak height of Lorentzian = 4 / (pi*hwhm)
% This defnition of the Lorentzian has Area = 1

area      = x(1);
hwhm      = x(2);
f0        = x(3);
phase     = x(4);
baseline0 = x(5);
baseline1 = x(6);

Absorption = 1/(2*pi) * area * hwhm ./ ((freq - f0).^2 + hwhm.^2) + ...
             1/(2*pi) * area * x(7) * hwhm ./ ((freq - f0 - 0.18).^2 + hwhm.^2);

Dispersion = 1/(2*pi) * area * (freq - f0) ./ ((freq - f0).^2 + hwhm.^2) + ...
             1/(2*pi) * area * x(7) * (freq - f0 - 0.18) ./ ((freq - f0 - 0.18).^2 + hwhm.^2);

F = cos(phase) * Absorption + sin(phase) * Dispersion + ...
    baseline0 + baseline1 * (freq - f0);
