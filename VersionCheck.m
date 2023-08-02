function varargout = VersionCheck(silent, currentVersion)
% Code adopted from Yair Altman's export_fig toolbox
% (https://github.com/altmany/export_fig)

% Check if there's a connection to the internet
try
    java.net.InetAddress.getByName('www.google.com');
catch
    warning('No internet connection. Skipping version check.');
    return
end

persistent lastCheckTime

if nargin < 2
    loadFile = which('GannetLoad');
    if isempty(loadFile)
        error('Cannot find GannetLoad.m; please ensure that your Gannet directory is included in your MATLAB search path.')
    end
    fileID = fopen(loadFile, 'rt');
    if fileID == -1
        error('Cannot read %s.', loadFile);
    end
    str = fread(fileID, Inf, '*uchar');
    fclose(fileID);
    str = char(str(:)');
    expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
    out = regexp(str, expression, 'names');
    currentVersion = out.version;
    if nargin < 1
        silent = 0;
    end
end

newVersionAvailable = 0;
if nargin < 2 || isempty(lastCheckTime) || (datetime('now') - lastCheckTime) > days(1)
    url = 'https://raw.githubusercontent.com/markmikkelsen/Gannet/main/GannetLoad.m';
    str = readURL(url);
    if isempty(str)
        if ~silent
            warning('Unable to check for Gannet updates; this may happen if your system has limited internet access.');
        end
        newVersionAvailable = 0;
    else
        expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
        out = regexp(str, expression, 'names');
        latestVersion = out.version;
        if str2double(latestVersion(regexpi(latestVersion,'\d'))) > str2double(currentVersion(regexpi(currentVersion,'\d')))
            newVersionAvailable = 1;
            msg = ['\n', ...
                   '***********************************************************************************************\n', ...
                   'A newer version of Gannet (%s) is available. You are currently using version %s.\n' ...
                   'You can download the newer version from GitHub or run UpdateGannet to install it directly.\n', ...
                   '***********************************************************************************************\n'];
            msg = hyperlink('https://github.com/markmikkelsen/Gannet', 'GitHub', msg);
            msg = hyperlink('matlab:UpdateGannet', 'UpdateGannet', msg);
            if ~silent
                fprintf(msg, latestVersion, currentVersion);
            end
        end
        lastCheckTime = datetime('now');
    end
end

if nargout == 1
    varargout{1} = currentVersion;
elseif nargout == 2
    varargout{1} = currentVersion;
    varargout{2} = newVersionAvailable;
end

    function str = readURL(url)
        % ARC 2023-06-27, exception handling also for urlread calls
        try
            if verLessThan('matlab','8.0')
                str = urlread(url); %#ok<URLRD> % R2012a or older (no timeout parameter)
            elseif verLessThan('matlab','8.4')
                str = urlread(url, 'Timeout', 5); %#ok<URLRD> % pre R2014b, no webread function
            else
                wo = weboptions('Timeout',5); % (note, 5 sec is consistent with Matlab default)
                str = char(webread(url, wo));
            end
            if size(str,1) > 1  % ensure a row-wise string
                str = str';
            end
        catch err
            if ~isempty(strfind(err.message,'404'))
                rethrow(err);
            else
                warning(err.message);
                str = '';
            end
        end
    end

end



