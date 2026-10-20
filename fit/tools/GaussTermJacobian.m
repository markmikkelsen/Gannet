function [G, dA, dw, dc] = GaussTermJacobian(A, w, c, freq)
% Gaussian term G = A*exp(w*(freq-c).^2) and its partial derivatives with
% respect to amplitude (A), width (w = -1/(2*sigma^2)) and center freq (c)

d  = freq - c;
dA = exp(w * d.^2);
G  = A * dA;
dw = G .* d.^2;
dc = -2 * w * G .* d;
