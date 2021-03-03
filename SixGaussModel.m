function F = SixGaussModel(x,freq)

% x(1) = Gaussian amplitude
% x(2) = 1/(2*sigma^2)
% x(3) = centre freq of peak
% x(4-6) = Gaussian 2
% x(7-9) = Gaussian 3
% x(10-12) = Gaussian 4
% x(13-15) = Gaussian 5
% x(16-18) = Gaussian 5
% x(19) = offset
% x(20) = slope
% x(21) = quadratic

F = x(1)*exp(x(2)*(freq-x(3)).*(freq-x(3))) + ...
    x(4)*exp(x(5)*(freq-x(6)).*(freq-x(6))) + ...
    x(7)*exp(x(8)*(freq-x(9)).*(freq-x(9))) + ...
    x(10)*exp(x(11)*(freq-x(12)).*(freq-x(12))) + ...
    x(13)*exp(x(14)*(freq-x(15)).*(freq-x(15))) + ...
    x(16)*exp(x(17)*(freq-x(18)).*(freq-x(18))) + ...
    x(19)+x(20)*(freq-x(3))+x(21)*(freq-x(3)).^2;