function A = SDATread(filename, da_xres)

% Open file to read reference scan data.
fid = fopen(filename, 'rb', 'ieee-le');
if fid == -1
    sprintf('Unable to locate file %s', filename);
    return
end

% Set up a structure to take the data:
totalframes = 1;
totalpoints = totalframes * da_xres * 2;

raw_data = freadVAXG(fid, totalpoints, 'float');

TempData = reshape(raw_data, [2 da_xres]);
A = squeeze(complex(TempData(1,:), TempData(2,:)));
fclose(fid);

end