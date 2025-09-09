function fd = dicom_open(strFilename)

%
% Gets a file descriptor to a DICOM file and check
% it vaguely for DICOM like properties...
%

fd = fopen(strFilename, 'rb');

if fd == -1
    fprintf('\nFailed to open DICOM file.');
    return;
end

% loose the header
fread(fd, 128, 'uchar');

if strcmp(sprintf('%c', fread(fd, 4, 'schar')),'DICM')
    %fprintf('\nFile appears to be valid DICOM.');
else
    fprintf('\nFile does NOT appear to be valid DICOM.');
    fclose(fd);
    fd = -1;
    return;
end
