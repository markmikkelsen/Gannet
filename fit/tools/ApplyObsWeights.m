function [F, J] = ApplyObsWeights(model, sqrtw, x, freq)
% Evaluates a signal model scaled by observation weights sqrt(w), and, if
% requested, the correspondingly weighted analytic Jacobian. Use as:
%   weightedModel = @(x,freq) ApplyObsWeights(@SomeModel, sqrt(w), x, freq);
% (an anonymous function forwards nargout, so the Jacobian is passed through)

if nargout > 1
    [F, J] = model(x, freq);
    J = sqrtw(:) .* J;
else
    F = model(x, freq);
end
F = reshape(sqrtw, size(F)) .* F;
