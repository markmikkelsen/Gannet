function [beta_hat, residual, h_tmp, resnorm, exitflag, output, lambda, jacobian] = ...
    FitSignalModel(model, freq, spec, baseline, beta0, lb, ub, lsqnlinopts, debug)

if nargin < 9
    debug = 0;
end

freq = freq(:);
spec = spec(:);
baseline = baseline(:);

% Function for problem solver
objFun = @(beta) SolveProblem(beta, freq, spec, baseline, model);

% Run nonlinear least squares optimization
[beta_hat, resnorm, residual, exitflag, output, lambda, jacobian] = ...
    lsqnonlin(objFun, beta0, lb, ub, lsqnlinopts);

if exitflag == -2 && ~debug
    error('Fitting failure! Set MRS_struct.p.debug = 1 in GannetPreInitialise.m for details.');
elseif exitflag == -2
    failure_dir = fullfile(pwd, 'Gannet_model_output', 'fail');
    if ~exist(failure_dir, 'dir')
        mkdir(failure_dir);
    end
    model_name_tok = regexp(func2str(model), '\w*Model\w*', 'match', 'once');
    if isempty(model_name_tok)
        model_name = 'UnknownModel';
    else
        model_name = model_name_tok;
    end
    failure_file = fullfile(failure_dir, ...
        sprintf('FitSignalModel_failure_%s_%s.mat', model_name, datetime('now', 'Format', 'yyMMdd_HHmmss')));
    model_str = func2str(model);
    bounds_check.lb_gt_ub    = find(lb > ub);
    bounds_check.beta0_lt_lb = find(beta0 < lb);
    bounds_check.beta0_gt_ub = find(beta0 > ub);
    call_stack = dbstack('-completenames');
    save(failure_file, ...
        'beta_hat', 'resnorm', 'residual', 'exitflag', 'output', 'lambda', 'jacobian', ...
        'model', 'model_str', 'freq', 'spec', 'baseline', 'beta0', 'lb', 'ub', 'lsqnlinopts', ...
        'bounds_check', 'call_stack');
    error(['Fitting failure! ' output.message ' lsqnonlin output saved to ' failure_file '.']);
end

h_tmp = figure('Visible', 'off');
% h_tmp = figure(333);
clf(h_tmp);
hold on;
plot(freq, spec, 'k', 'LineWidth', 1);
plot(freq, model(beta_hat, freq) + baseline, 'r', 'LineWidth', 1);
plot(freq, baseline, 'LineWidth', 1, 'Color', '#FCAF0A');
plot(freq, residual - 0.25, 'k');
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
