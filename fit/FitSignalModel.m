function [beta_hat, residual, h_tmp, resnorm, exitflag, output, lambda, jacobian] = ...
    FitSignalModel(model, freq, data, baseline, beta0, lb, ub, lsqnlinopts)

freq = freq(:);
data = data(:);
baseline = baseline(:);

% Function for problem solver
objFun = @(beta) SolveProblem(beta, freq, data, baseline, model);

% Run nonlinear least squares optimization
[beta_hat, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objFun, beta0, lb, ub, lsqnlinopts);

h_tmp = figure('Visible', 'off');
% h_tmp = figure(333);
clf(h_tmp);
hold on;
plot(freq, data, 'k', 'LineWidth', 1);
plot(freq, model(beta_hat, freq) + baseline, 'r', 'LineWidth', 1);
plot(freq, baseline);
plot(freq, residual, 'k');
hold off;
xlabel('ppm','FontSize',16);
set(gca,'XDir','reverse','TickDir','out');
legend({'data','model + baseline','baseline','residual'}, ...
    'Box','off','Location','best');
drawnow;

end


function r = SolveProblem(beta, freq, data, baseline, model)

% 1) Data fit residuals
y_hat = model(beta, freq) + baseline;
r = data(:) - y_hat(:);

% 2) Parameter constraint term based on baseline
% Weight by sqrt(lambda) so lambda acts like a penalty weight
% r_constraint = sqrt(lambda) * (y_hat - baseline);
% r_constraint = sqrt(lambda) * baseline;

% Stack into a single residual vector
% r = [r_data; r_constraint];
% r = r_data;

end
