function F = GABAGlxModel(x,freq)
% Function for GABA+Glx model

% Three Gaussians
%  x(1)  = gaussian amplitude 1
%  x(2)  = width 1 ( 1/(2*sigma^2) )
%  x(3)  = centre freq of peak 1
%  x(4)  = gaussian amplitude 2
%  x(5)  = width 2 ( 1/(2*sigma^2) )
%  x(6)  = centre freq of peak 2
%  x(7)  = gaussian amplitude 3
%  x(8)  = width 3 ( 1/(2*sigma^2) )
%  x(9)  = centre freq of peak 3
%  x(10) = linear baseline slope
%  x(11) = sine baseline term
%  x(12) = cosine baseline term

% MM: Allowing peaks to vary individually seems to work better than keeping
% the distance fixed (i.e., including J in the function)

F = x(1) * exp(x(2) * (freq - x(3)).^2) + ...
    x(4) * exp(x(5) * (freq - x(6)).^2) + ...
    x(7) * exp(x(8) * (freq - x(9)).^2) + ...
    x(10) * (freq - x(3)) + ...
    x(11) * sin(pi * freq/1.31/4) + x(12) * cos(pi * freq/1.31/4);
