function F = ThreeLorentzModel_phased_noBaseline(x,freq)
% ThreeLorentzModel with phase and no baseline
% Based on ﻿Marshall & Roe, Anal Chem, 1978;50(6):756-763,
% doi:10.1021/ac50027a023

H   = x(1); % amplitude (center peak)
a   = x(2); % amplitude scaling factor
b   = x(3); % amplitude scaling factor
T2  = x(4); % T2 relaxation time constant (in ms)
f0  = x(5); % frequency (in ppm)
J   = x(6); % J-coupling constant (in ppm)
phi = x(7); % phase (in rad)

A1 = cos(phi) .* ((a .* H .* T2) ./ (1 + ((f0 + J) - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((a .* H .* ((f0 + J) - freq) .* T2.^2) ./ (1 + ((f0 + J) - freq).^2 .* T2.^2));

A2 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

A3 = cos(phi) .* ((b .* H .* T2) ./ (1 + ((f0 - J) - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((b .* H .* ((f0 - J) - freq) .* T2.^2) ./ (1 + ((f0 - J) - freq).^2 .* T2.^2));

F = A1 + A2 + A3;
