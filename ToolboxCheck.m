function ToolboxCheck

req_toolboxes = {'Image Processing Toolbox'
                 'Optimization Toolbox'
                 'Signal Processing Toolbox'
                 'Statistics and Machine Learning Toolbox'};

instld_toolboxes = ver;
instld_toolboxes = {instld_toolboxes.Name}';
missing          = [];

for ii = 1:length(req_toolboxes)
    if ~any(strcmp(instld_toolboxes, req_toolboxes(ii)))
        missing = [missing; req_toolboxes(ii)]; %#ok<*AGROW>
    end
end

if ~isempty(missing)
    msg = 'The following MATLAB toolboxes are missing and need to be installed to run Gannet:\n\n';
    for ii = 1:length(missing)
        msg = [msg missing{ii} '\n'];
    end
    msg = [msg '\n'];
    error(sprintf(msg));
end



