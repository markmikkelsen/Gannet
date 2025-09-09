function UpdateGannet
% Code adopted from Yair Altman's export_fig toolbox
% (https://github.com/altmany/export_fig) and SPM12's update tool

% Check if there's a connection to the internet
try
    java.net.InetAddress.getByName('www.google.com');
catch
    error('No internet connection. Cannot run UpdateGannet.');
end

% First, check if a new version of Gannet is available; exit otherwise
[currentVersion, newVersionAvailable] = VersionCheck(1);
if ~newVersionAvailable
    fprintf('\nYour version of Gannet (%s) is the latest version.\n\n', currentVersion);
    return
end

% Present a warning to the user and ask if update should proceed
opts.Default = 'No';
opts.Interpreter = 'tex';
answer = questdlg({['\fontsize{12}WARNING: Running the automated updater will permanently ' ...
                    'replace the contents of the current Gannet folder in your search ' ...
                    'path with the latest code from GitHub.'] ...
                    '' ...
                   ['{\color{red}Do NOT use this updater if you cloned the Gannet GitHub repository ' ...
                    'using git or GitHub Desktop. Pull the latest commits in the usual way instead.}'] ...
                    '' ...
                    'Do you wish to continue?'}, 'Update Gannet', 'Yes', 'No', opts);
switch answer
    case 'Yes'
        fprintf('\nUpdating Gannet...\n\n');
    case 'No'
        fprintf('\nExiting updater...\n\n');
        return
end

% Remove the current Gannet folder from the search path
gannetPath = fileparts(which(mfilename('fullpath')));
searchPath = textscan(path, '%s', 'delimiter', pathsep);
searchPath = searchPath{1};
i = strncmp(gannetPath, searchPath, length(gannetPath));
searchPath(i) = [];
searchPath = strcat(searchPath, pathsep);
path(strcat(searchPath{:}));

% Download the latest version of Gannet
zipURL = 'https://github.com/markmikkelsen/Gannet/archive/refs/heads/main.zip';
targetFolder = fullfile(pwd, ['tmp_' randsample(['A':'Z','0':'9'],5)]);
mkdir(targetFolder);
targetFilename = fullfile(targetFolder, [char(datetime('now','Format','y-MM-dd')) '.zip']);
websave(targetFilename, zipURL);
newFilenames = unzip(targetFilename, targetFolder);

% Delete the old Gannet folder in the search path and replace it with the new one
rmdir(gannetPath,'s');
movefile(newFilenames{1}, gannetPath);
rmdir(targetFolder,'s');

% Add the new Gannet folder and its subfolders to search path
addpath(genpath(gannetPath));

% Notify the user and rehash
url = 'https://raw.githubusercontent.com/markmikkelsen/Gannet/main/GannetLoad.m';
str = readURL(url);
expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
latestVersion = out.version;
fprintf('\nSuccessfully updated Gannet to version %s!\n\n', latestVersion);
rehash;

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


