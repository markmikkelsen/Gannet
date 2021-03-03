function CheckTargets(MRS_struct)

filepath = fullfile(fileparts(which(mfilename('fullpath'))), 'GannetPreInitialise.m');
msg      = 'Incorrect targets entered in GannetPreInitialise.m. Check spelling, allowable options and order of targets.';
msg      = hyperlink(['matlab: opentoline(''' filepath ''', 4, 0)'], 'Incorrect targets entered in GannetPreInitialise.m', msg);

switch num2str(length(MRS_struct.p.target))
    case '1'
        if any(strcmp(MRS_struct.p.target,{'GABA','Glx','GABAGlx','GSH','Lac','EtOH'}))
            if MRS_struct.p.phantom && strcmp(MRS_struct.p.target,'GABAGlx')
                error(msg);
            end
            if MRS_struct.p.HERMES
                msg = 'MRS_struct.p.HERMES is set to 1 in GannetPreInitialise.m. Add a second target metabolite or set flag to 0.';
                msg = hyperlink(['matlab: opentoline(''' filepath ''', 4, 0)'], 'Add a second target metabolite', msg);
                msg = hyperlink(['matlab: opentoline(''' filepath ''', 34, 0)'], 'set flag to 0', msg);
                error(msg);
            end
        else
            PassError;
        end
    case '2'
        if any([all(strcmp(MRS_struct.p.target,{'GABAGlx','GSH'})) ...
                all(strcmp(MRS_struct.p.target,{'GABA','GSH'})) ...
                all(strcmp(MRS_struct.p.target,{'Glx','GSH'})) ...
                all(strcmp(MRS_struct.p.target,{'Lac','GSH'}))])
            if MRS_struct.p.phantom && any(strcmp(MRS_struct.p.target,'GABAGlx'))
                error(msg);
            end
            if ~MRS_struct.p.HERMES
                msg = 'Two target metabolites detected. MRS_struct.p.HERMES must be set to 1 in GannetPreInitialise.m.';
                msg = hyperlink(['matlab: opentoline(''' filepath ''', 34, 0)'], 'MRS_struct.p.HERMES must be set to 1', msg);
                error(msg);
            end
        else
            error(msg);
        end
    case '3'
        if all(strcmp(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            if ~MRS_struct.p.HERMES
                msg = 'Three target metabolites detected. MRS_struct.p.HERMES must be set to 1 in GannetPreInitialise.m.';
                msg = hyperlink(['matlab: opentoline(''' filepath ''', 34, 0)'], 'MRS_struct.p.HERMES must be set to 1', msg);
                error(msg);
            end
        else
            error(msg);
        end
end

end



