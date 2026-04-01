function F = EtOHModel_noBaseline(x, freq)
% Function for EtOH model with no baseline

F = x(1) ./ (1 + ((freq - x(2)) / (x(3)/2)).^2) + ...
    x(4) ./ (1 + ((freq - x(5)) / (x(6)/2)).^2);
