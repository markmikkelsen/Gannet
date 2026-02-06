function F = ThreeLorentzModel_phased_noBaseline(x,freq)
% ThreeLorentzModel with phase and no baseline
% Based on ﻿Marshall & Roe, Anal Chem, 1978;50(6):756-763,
% doi:10.1021/ac50027a023

H   = x(1);  % amplitude
T2  = x(2);  % T2 relaxation time constant
f0  = x(3);  % frequency (in ppm)
phi = x(4);  % phase (in rad)

A1 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

A2 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

A3 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

F = A1 + A2 + A3;
