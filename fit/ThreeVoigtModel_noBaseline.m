function F = ThreeVoigtModel_noBaseline(x, freq)

% Formula taken from Marshall et al. Magn Reson Med. 1997;37(5):651-657
% doi:10.1002/mrm.1910370504

r1 = x(1);
a1 = x(2);
f1 = x(3);
W1 = x(4);

L1 = a1 ./ (1 + ((freq - f1) ./ (W1./2)).^2);
G1 = a1 .* exp(-log(2) .* ((freq - f1) ./ (W1./2)).^2);
V1 = r1 .* L1 + (1 - r1) .* G1;

r2 = x(5);
a2 = x(6);
f2 = x(7);
W2 = x(8);

L2 = a2 ./ (1 + ((freq - f2) ./ (W2./2)).^2);
G2 = a2 .* exp(-log(2) .* ((freq - f2) ./ (W2./2)).^2);
V2 = r2 .* L2 + (1 - r2) .* G2;

r3 = x(9);
a3 = x(10);
f3 = x(11);
W3 = x(12);

L3 = a3 ./ (1 + ((freq - f3) ./ (W3./2)).^2);
G3 = a3 .* exp(-log(2) .* ((freq - f3) ./ (W3./2)).^2);
V3 = r3 .* L3 + (1 - r3) .* G3;

F = V1 + V2 + V3;
