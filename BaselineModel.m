function F = BaselineModel(x,freq)
% Function for Baseline Model

F = x(2) * (freq - x(1)) + x(3);
