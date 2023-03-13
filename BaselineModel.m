function F = BaselineModel(x, freq)
% Function for baseline model

F = x(2) * (freq - x(1)) + x(3);
