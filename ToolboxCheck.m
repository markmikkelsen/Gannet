function ToolboxCheck

required_toolboxes = {'Image Processing Toolbox'
                      'Optimization Toolbox'
                      'Signal Processing Toolbox'
                      'Statistics and Machine Learning Toolbox'};

installed_toolboxes = ver;
installed_toolboxes = {installed_toolboxes.Name}';
missing             = [];

for ii = 1:length(required_toolboxes)
    if ~any(strcmp(installed_toolboxes, required_toolboxes(ii)))
        missing = [missing; required_toolboxes(ii)]; %#ok<*AGROW> 
    end
end

if ~isempty(missing)
    msg = 'The following MATLAB toolboxes are missing and need to be installed to run Gannet (see instructions):\n\n';
    for ii = 1:length(missing)
        msg = [msg missing{ii} '\n'];
    end
    msg(end-1:end) = [];
    msg = hyperlink('https://www.mathworks.com/matlabcentral/answers/101885-how-do-i-install-additional-toolboxes-into-an-existing-installation-of-matlab', ...
                    'see instructions', msg);
    fprintf('\n');
    error(sprintf(msg));
end



