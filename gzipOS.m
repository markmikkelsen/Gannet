function gzipOS(fname)
% Use pigz or system gzip if available (faster)

persistent cmd; % command to gzip
if isempty(cmd)
    cmd = check_gzip('gzip');
    if ischar(cmd)
        cmd = @(nam){cmd '-nf' nam};
    elseif islogical(cmd) && ~cmd
        fprintf(2, ['None of system pigz, gzip or Matlab gzip available. ' ...
            'Files are not compressed into gz.\n']);
    end
end

if islogical(cmd)
    if cmd, gzip(fname); deleteFile(fname); end
    return;
end

[err, str] = jsystem(cmd(fname));
if err && ~exist(strcat(fname, '.gz'), 'file')
    try
        gzip(fname); deleteFile(fname);
    catch
        fprintf(2, 'Error during compression: %s\n', str);
    end
end


function cmd = check_gzip(gz_unzip)
% Deal with pigz/gzip on path or in nii_tool folder, and matlab gzip/gunzip

m_dir = fileparts(mfilename('fullpath'));
% away from pwd, so use OS pigz if both exist. Avoid error if pwd changed later
if strcmpi(pwd, m_dir), cd ..; clnObj = onCleanup(@() cd(m_dir)); end
if isunix
    pth1 = getenv('PATH');
    if isempty(strfind(pth1, '/usr/local/bin')) %#ok<*STREMP>
        pth1 = [pth1 ':/usr/local/bin'];
        setenv('PATH', pth1);
    end
end

% first, try system pigz
[err, ~] = jsystem({'pigz' '-V'});
if ~err, cmd = 'pigz'; return; end

% next, try pigz included with nii_tool
cmd = [m_dir '/pigz'];
if ismac % pigz for mac is not included in the package
    if strcmp(gz_unzip, 'gzip')
        fprintf(2, [' Please install pigz for fast compression: ' ...
            'http://macappstore.org/pigz/\n']);
    end
elseif isunix % linux
    [st, val] = fileattrib(cmd);
    if st && ~val.UserExecute, fileattrib(cmd, '+x'); end
end

[err, ~] = jsystem({cmd '-V'});
if ~err, return; end

% Third, try system gzip/gunzip
[err, ~] = jsystem({gz_unzip '-V'}); % gzip/gunzip on system path?
if ~err, cmd = gz_unzip; return; end

% Lastly, use Matlab gzip/gunzip if java avail
cmd = usejava('jvm');


function [err, out] = jsystem(cmd)
% faster than system: based on https://github.com/avivrosenberg/matlab-jsystem

% cmd is cell str, no quotation marks needed for file names with space.
cmd = cellstr(cmd);
try
    pb = java.lang.ProcessBuilder(cmd);
    pb.redirectErrorStream(true); % ErrorStream to InputStream
    process = pb.start();
    scanner = java.util.Scanner(process.getInputStream).useDelimiter('\A');
    if scanner.hasNext(), out = char(scanner.next()); else, out = ''; end
    err = process.exitValue; % err = process.waitFor() may hang
    if err, error('java.lang.ProcessBuilder error'); end
catch % fallback to system() if java fails like for Octave
    cmd = regexprep(cmd, '.+? .+', '"$0"'); % double quotes if with middle space
    [err, out] = system(sprintf('%s ', cmd{:}, '2>&1'));
end


function deleteFile(fname)
% Delete file in background

if ispc, system(['start "" /B del "' fname '"']);
else, system(['rm "' fname '" &']);
end



