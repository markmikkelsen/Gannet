function F = ThreeGaussModel_noBaseline(x, freq)
% Function for three-Gaussian model with no baseline

% x(1-3) = amplitudes
% x(4-6) = widths (1/(2*sigma^2))
% x(7-9) = center freqs

F = x(1) * exp(x(4) * (freq - x(7)).^2) + ...
    x(2) * exp(x(5) * (freq - x(8)).^2) + ...
    x(3) * exp(x(6) * (freq - x(9)).^2);
