function [F, J] = ThreeLorentzModel_phased_noBaseline(x,freq)
% ThreeLorentzModel with phase and no baseline
% Based on Marshall & Roe, Anal Chem, 1978;50(6):756-763,
% doi:10.1021/ac50027a023
%
% If requested, J returns the analytic Jacobian dF/dx (numel(freq) x 7)

H   = x(1); % amplitude (center peak)
a   = x(2); % amplitude scaling factor
b   = x(3); % amplitude scaling factor
T2  = x(4); % T2 relaxation time constant (in ms)
f0  = x(5); % frequency (in ppm)
J   = x(6); % J-coupling constant (in ppm)
phi = x(7); % phase (in rad)

A1 = cos(phi) .* ((a .* H .* T2) ./ (1 + ((f0 + J) - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((a .* H .* ((f0 + J) - freq) .* T2.^2) ./ (1 + ((f0 + J) - freq).^2 .* T2.^2));

A2 = cos(phi) .* ((H .* T2) ./ (1 + (f0 - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((H .* (f0 - freq) .* T2.^2) ./ (1 + (f0 - freq).^2 .* T2.^2));

A3 = cos(phi) .* ((b .* H .* T2) ./ (1 + ((f0 - J) - freq).^2 .* T2.^2))  - ...
    sin(phi) .* ((b .* H .* ((f0 - J) - freq) .* T2.^2) ./ (1 + ((f0 - J) - freq).^2 .* T2.^2));

F = A1 + A2 + A3;

if nargout > 1
    % Each peak: Ak = sk*H*N/D, with N = cos(phi)*T2 - sin(phi)*d*T2^2,
    % D = 1 + d^2*T2^2, d = ck - freq
    f  = freq(:);
    cp = cos(phi);
    sp = sin(phi);
    scale  = [a 1 b];
    center = [f0+J f0 f0-J];
    dcdJ   = [1 0 -1]; % d(center)/dJ
    Jac = zeros(numel(f), numel(x));
    for k = 1:3
        d   = center(k) - f;
        D   = 1 + d.^2 * T2^2;
        N   = cp*T2 - sp*d*T2^2;
        sH  = scale(k) * H;
        dAdd = sH * (-sp*T2^2 .* D - N .* 2.*d*T2^2) ./ D.^2; % dAk/dd
        Jac(:,1) = Jac(:,1) + scale(k) * N ./ D;
        Jac(:,4) = Jac(:,4) + sH * ((cp - 2*sp*d*T2) .* D - N .* 2.*d.^2*T2) ./ D.^2;
        Jac(:,5) = Jac(:,5) + dAdd;
        Jac(:,6) = Jac(:,6) + dcdJ(k) * dAdd;
        Jac(:,7) = Jac(:,7) + sH * (-sp*T2 - cp*d*T2^2) ./ D;
    end
    Jac(:,2) = H * (cp*T2 - sp*(center(1)-f)*T2^2) ./ (1 + (center(1)-f).^2 * T2^2);
    Jac(:,3) = H * (cp*T2 - sp*(center(3)-f)*T2^2) ./ (1 + (center(3)-f).^2 * T2^2);
    J = Jac;
end
