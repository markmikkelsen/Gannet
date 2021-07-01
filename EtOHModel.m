function F = EtOHModel(x,freq)
% Function for EtOH model

L1 = x(1) ./ (1 + ((freq - x(2)) / (x(3)/2)).^2);
L2 = x(4) ./ (1 + ((freq - x(5)) / (x(6)/2)).^2);
B  = x(7) .* (freq - x(3)) + x(8);
F  = L1 + L2 + B;