function [v, w] = WACFM(x, costFun)
% Weighted averaging based on criterion function minimization (WACFM).
% Algorithm from Pander T. A new approach to robust, weighted signal
% averaging. Biocybern Biomed Eng. 2015;35(4):317-327.
% doi:10.1016/j.bbe.2015.06.002

[~,N]  = size(x);
kStop  = 100;
v      = median(x,2);
w      = zeros(kStop,N);
e      = 1e-6;
m      = 2;
p      = 0.2;
sigma  = 1;
const  = max(abs(v)); % constant to increase robustness to local minima
                      % (﻿Kotowski et al. Biocybern Biomed Eng. 2019.
                      % doi:10.1016/j.bbe.2019.09.002)

for k = 1:kStop
    z = x - v;
    w(k,:) = weights(z, costFun, const);

    % figure(23);
    % cla;
    % plot(w(1:k,:)');
    % drawnow;
    % pause(0.25);

    if (k > 1 && norm(w(k,:) - w(k-1,:)) < e) || k == kStop
        w = w(k,:);
        switch costFun % normalize optimal weights so they sum to unity
            case 'square'
                w = w.^m ./ sum(w.^m);
            case 'GCD'
                w = w.^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1 ./ ...
                    sum(w.^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1);
        end
        break
    end

    switch costFun
        case 'square'
            v = sum(w(k,:).^m .* x,2) ./ sum(w(k,:).^m);
        case 'GCD'
            v = sum(w(k,:).^m .* x .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1,2) ./ ...
                sum(w(k,:).^m .* sum((abs(z).^(2-p)) ./ p .* (sigma.^p + abs(z).^p)).^-1);
    end
end

    function w = weights(z, costFun, const)
        switch costFun
            case 'square'
                w = (vecnorm(z) + const).^(2./(1-m)) ./ ...
                    sum((vecnorm(z) + const).^(2./(1-m)));
            case 'GCD'
                w = sum(log(1 + (abs(z) ./ sigma).^p) + const).^(1./(1-m)) ./ ...
                    sum(sum(log(1 + (abs(z) ./ sigma).^p) + const).^(1./(1-m)));
        end
    end

end
