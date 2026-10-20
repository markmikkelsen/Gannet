function V = VoigtModel(x,freq)
% Formula taken from Marshall et al. Magn Reson Med. 1997;37(5):651-657
% doi:10.1002/mrm.1910370504

r  = x(1);
a  = x(2);
f0 = x(3);
W  = x(4);
m  = x(5);
b  = x(6);

L = a ./ (1 + ((freq - f0) ./ (W./2)).^2);
G = a .* exp(-log(2) .* ((freq - f0) ./ (W./2)).^2);
V = exp(-abs(r)) .* L + (1 - exp(-abs(r))) .* G + m .* (freq - f0) + b;
