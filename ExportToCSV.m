function MRS_struct = ExportToCSV(MRS_struct, vox, module)

round2 = @(x) round(x*1e3)/1e3;

if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
    n_rep = [size(MRS_struct.metabfile,2)/2 1];
else
    n_rep = [size(MRS_struct.metabfile,2) 1];
end
out.MATLAB_ver       = cellstr(repmat(version('-release'), n_rep));
out.Gannet_ver       = cellstr(repmat(MRS_struct.version.Gannet, n_rep));
out.date_of_analysis = cellstr(repmat(char(datetime('now','Format','y-MM-dd')), n_rep));


%%% 1. Extract data from GannetFit %%%

if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
    filename = MRS_struct.metabfile(1,1:2:end)';
else
    filename = MRS_struct.metabfile(1,:)';
end
for ii = 1:length(filename)
    [~,b,c] = fileparts(filename{ii});
    out.filename(ii,1) = cellstr([b c]);
end
out.avg_delta_F0 = MRS_struct.out.AvgDeltaF0(:);

metabs = {'GABA','Glx','GSH','EtOH','Lac','water','Cr','Cho','NAA'};

for ii = 1:length(metabs)
    if ~isfield(MRS_struct.out.(vox), metabs{ii})
        continue
    end
    out.(metabs{ii}).area      = MRS_struct.out.(vox).(metabs{ii}).Area(:);
    out.(metabs{ii}).FWHM      = MRS_struct.out.(vox).(metabs{ii}).FWHM(:);
    out.(metabs{ii}).SNR       = MRS_struct.out.(vox).(metabs{ii}).SNR(:);
    out.(metabs{ii}).fit_error = MRS_struct.out.(vox).(metabs{ii}).FitError(:);
    if ~strcmp(metabs{ii}, 'water')
        if ~strcmp(metabs{ii}, 'Cr')
            out.(metabs{ii}).fit_error_Cr = MRS_struct.out.(vox).(metabs{ii}).FitError_Cr(:);
        end
        if strcmp(MRS_struct.p.reference, 'H2O')
            out.(metabs{ii}).fit_error_w = MRS_struct.out.(vox).(metabs{ii}).FitError_W(:);
        end
        if ~strcmp(metabs{ii}, 'Cr')
            out.(metabs{ii}).conc_Cr = MRS_struct.out.(vox).(metabs{ii}).ConcCr(:);
        end
        if strcmp(MRS_struct.p.reference, 'H2O')
            out.(metabs{ii}).conc_iu = MRS_struct.out.(vox).(metabs{ii}).ConcIU(:);
        end
    end
end

T = table(out.MATLAB_ver, out.Gannet_ver, out.date_of_analysis, out.filename, round2(out.avg_delta_F0), ...
    'VariableNames', {'MATLAB_version', 'Gannet_version', 'date_of_analysis', 'filename', 'avg_delta_F0'});

field_names = fieldnames(out);

for ii = 1:length(field_names)
    if any(strcmp(field_names{ii}, metabs))
        sub_field_names = fieldnames(out.(field_names{ii}));
        for jj = 1:length(sub_field_names)
            if strcmp(sub_field_names{jj}, 'area')
                U = table(out.(field_names{ii}).(sub_field_names{jj}), ...
                    'VariableNames', {[field_names{ii} '_' sub_field_names{jj}]});
            else
                U = table(round2(out.(field_names{ii}).(sub_field_names{jj})), ...
                    'VariableNames', {[field_names{ii} '_' sub_field_names{jj}]});
            end
            T = [T U]; %#ok<*AGROW>
        end
    end
end

% Create CSV filename
if isfield(MRS_struct.out.(vox), 'csv_name')
    csv_name = MRS_struct.out.(vox).csv_name;
else
    csv_name = fullfile(pwd, 'Gannet_output.csv');
    if exist(csv_name, 'file')
        run_count = 1;
        csv_name = fullfile(pwd, ['Gannet_output' num2str(run_count) '.csv']);
        while 1
            if exist(csv_name, 'file')
                run_count = run_count + 1;
                csv_name  = fullfile(pwd, ['Gannet_output' num2str(run_count) '.csv']);
            else
                break
            end
        end
    end
    MRS_struct.out.(vox).csv_name = csv_name;
end

% End if function invoked in GannetFit
if strcmp(module, 'fit')
    % Convert empty cells into NaNs
    for ii = 1:size(T,2)
        if ~iscell(T(:,ii).(T.Properties.VariableNames{ii}))
            T(~T(:,ii).(T.Properties.VariableNames{ii}),ii) = {NaN};
        end
    end
    fprintf('\nExporting results to %s\n', [csv_name '...']);
    writetable(T, csv_name);
    return
end


%%% 2. Extract data from GannetSegment %%%

out.tissue.fGM  = MRS_struct.out.(vox).tissue.fGM(:);
out.tissue.fWM  = MRS_struct.out.(vox).tissue.fWM(:);
out.tissue.fCSF = MRS_struct.out.(vox).tissue.fCSF(:);

metabs = {'GABA','Glx','GSH','EtOH','Lac','Cr','Cho','NAA'};

if strcmp(MRS_struct.p.reference, 'H2O')
    for ii = 1:length(metabs)
        if ~isfield(MRS_struct.out.(vox), metabs{ii})
            continue
        end
        out.(metabs{ii}).ConcIU_CSFcorr = MRS_struct.out.(vox).(metabs{ii}).ConcIU_CSFcorr(:);
    end
end

field_names = fieldnames(out);
X = table;
V = table;

for ii = 1:length(field_names)
    if any(strcmp(field_names{ii}, metabs)) && strcmp(MRS_struct.p.reference, 'H2O')
        sub_field_names = fieldnames(out.(field_names{ii}));
        Y = table(round2(out.(field_names{ii}).(sub_field_names{end})), ...
            'VariableNames', {[field_names{ii} '_' sub_field_names{end}]});
        X = [X Y];
    elseif strcmp(field_names{ii}, 'tissue')
        sub_field_names = fieldnames(out.(field_names{ii}));
        for jj = 1:3
            U = table(round2(out.(field_names{ii}).(sub_field_names{jj})), ...
                'VariableNames', sub_field_names(jj));
            V = [V U];
        end
    end
end

T = [T V X]; % doing it this way so that tissue fractions come before the CSF-corrected values in the .csv file

% End if function invoked in GannetSegment
if strcmp(module, 'segment')
    % Convert empty cells into NaNs
    for ii = 1:size(T,2)
        if ~iscell(T(:,ii).(T.Properties.VariableNames{ii}))
            T(~T(:,ii).(T.Properties.VariableNames{ii}),ii) = {NaN};
        end
    end
    fprintf('\nUpdating results in %s\n', [csv_name '...']);
    writetable(T, csv_name);
    return
end


%%% 3. Extract data from GannetQuantify %%%

if strcmp(MRS_struct.p.reference, 'H2O')
    for ii = 1:length(metabs)
        if ~isfield(MRS_struct.out.(vox), metabs{ii})
            continue
        end
        out.(metabs{ii}).ConcIU_TissCorr              = MRS_struct.out.(vox).(metabs{ii}).ConcIU_TissCorr(:);
        out.(metabs{ii}).ConcIU_AlphaTissCorr         = MRS_struct.out.(vox).(metabs{ii}).ConcIU_AlphaTissCorr(:);
        out.(metabs{ii}).ConcIU_AlphaTissCorr_GrpNorm = MRS_struct.out.(vox).(metabs{ii}).ConcIU_AlphaTissCorr_GrpNorm(:);
        out.(metabs{ii}).alpha                        = MRS_struct.out.(vox).(metabs{ii}).alpha(:);
    end
end

field_names = fieldnames(out);

for ii = 1:length(field_names)
    if any(strcmp(field_names{ii}, metabs))
        sub_field_names = fieldnames(out.(field_names{ii}));
        for jj = length(sub_field_names)-3:length(sub_field_names)
            U = table(round2(out.(field_names{ii}).(sub_field_names{jj})), ...
                'VariableNames', {[field_names{ii} '_' sub_field_names{jj}]});
            T = [T U];
        end
    end
end

% Convert empty cells into NaNs
for ii = 1:size(T,2)
    if ~iscell(T(:,ii).(T.Properties.VariableNames{ii}))
        T(~T(:,ii).(T.Properties.VariableNames{ii}),ii) = {NaN};
    end
end
fprintf('\nUpdating results in %s\n', [csv_name '...']);
writetable(T, csv_name);



