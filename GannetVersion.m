function GannetVersion
% Return the version number of Gannet

loadFile = which('GannetLoad');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    error('Can''t read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nYour version of Gannet is %s\n\n', out.version);