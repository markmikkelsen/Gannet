function F = ThreeLorentzModel_noBaseline(x,freq)
% ThreeLorentzModel with phase and no baseline
% Based on ﻿Marshall & Roe, Anal Chem, 1978;50(6):756-763,
% doi:10.1021/ac50027a023

H   = x(1);  % amplitude of middle peak
a   = x(2);  % amplitude scaling factor for left outer peak
b   = x(3);  % amplitude scaling factor for right outer peak
T2  = x(4);  % T2 relaxation time constant
f0  = x(5);  % frequency of middle peak (in ppm)
J   = x(6);  % J-coupling constant (in ppm)
phi = x(7);  % phase (in deg)

A0 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
        sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

Aa = cos(phi) .* ((H .* a .* T2) ./ (1 + (f0 + J - freq).^2 .* T2.^2)) - ...
        sin(phi) .* ((H .* a .* (f0 + J - freq) .* T2.^2) ./ (1 + (f0 + J - freq).^2 .* T2.^2));

Ab = cos(phi) .* ((H .* b .* T2) ./ (1 + (f0 - J - freq).^2 .* T2.^2)) - ...
        sin(phi) .* ((H .* b .* (f0 - J - freq) .* T2.^2) ./ (1 + (f0 - J - freq).^2 .* T2.^2));

F = A0 + Aa + Ab;
