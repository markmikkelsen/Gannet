function GannetVersion
% Return the version numbers of Gannet and its modules

% GannetLoad
loadFile = which('GannetLoad');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    fclose(fileID);
    error('Cannot read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\n-------------------------------');
fprintf('\nYour version of Gannet is %s', out.version);
fprintf('\n-------------------------------');

expression = '(?<field>MRS_struct.version.load   = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nGannetLoad:          %s', out.version);

% GannetFit 
loadFile = which('GannetFit');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    fclose(fileID);
    error('Cannot read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.fit = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nGannetFit:           %s', out.version);

% GannetCoRegister 
loadFile = which('GannetCoRegister');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    fclose(fileID);
    error('Cannot read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.coreg = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nGannetCoRegister:    %s', out.version);

% GannetSegment 
loadFile = which('GannetSegment');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    fclose(fileID);
    error('Cannot read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.segment = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nGannetSegment:       %s', out.version);


% GannetQuantify 
loadFile = which('GannetQuantify');
fileID = fopen(loadFile, 'rt');
if fileID == -1
    fclose(fileID);
    error('Cannot read %s.', loadFile);
end
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.quantify = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
fprintf('\nGannetQuantify:      %s\n\n', out.version);



