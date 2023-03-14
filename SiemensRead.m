function MRS_struct = SiemensRead(MRS_struct, off_filename, on_filename, water_filename)

ii = MRS_struct.ii;

% First load in the OFF data
rda_filename = off_filename;
fid = fopen(rda_filename);

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';

tline = fgets(fid);

while isempty(strfind(tline, head_end_text))

    tline = fgets(fid);

    if isempty(strfind(tline, head_start_text)) + isempty(strfind(tline, head_end_text)) == 2

        % Store this data in the appropriate format
        occurence_of_colon = strfind(tline,':');
        variable = tline(1:occurence_of_colon-1);
        value    = tline(occurence_of_colon+1:length(tline));

        switch variable
            case {'PatientID', 'PatientName', 'StudyDescription', 'PatientBirthDate', 'StudyDate', 'StudyTime', 'PatientAge', 'SeriesDate', ...
                  'SeriesTime', 'SeriesDescription', 'ProtocolName', 'PatientPosition', 'ModelName', 'StationName', 'InstitutionName', ...
                  'DeviceSerialNumber', 'InstanceDate', 'InstanceTime', 'InstanceComments', 'SequenceName', 'SequenceDescription', 'Nucleus', ...
                  'TransmitCoil'}
                rda.(variable) = value;
            case {'PatientSex'}
                % Sex converter! (int to M,F,U)
                switch value
                    case 0
                        rda.sex = 'Unknown';
                    case 1
                        rda.sex = 'Male';
                    case 2
                        rda.sex = 'Female';
                end
            case {'SeriesNumber', 'InstanceNumber', 'AcquisitionNumber', 'NumOfPhaseEncodingSteps', 'NumberOfRows', 'NumberOfColumns', 'VectorSize'}
                %Integers
                rda.(variable) = str2double(value);
            case {'PatientWeight', 'TR', 'TE', 'TM', 'DwellTime', 'NumberOfAverages', 'MRFrequency', 'MagneticFieldStrength', 'FlipAngle', ...
                    'SliceThickness',  'FoVHeight', 'FoVWidth', 'PercentOfRectFoV', 'PixelSpacingRow', 'PixelSpacingCol'}
                %Floats
                rda.(variable) = str2double(value);
            case {'SoftwareVersion[0]'}
                rda.software_version = value;
            case {'CSIMatrixSize[0]'}
                rda.CSIMatrix_Size(1) = str2double(value);
            case {'CSIMatrixSize[1]'}
                rda.CSIMatrix_Size(2) = str2double(value);
            case {'CSIMatrixSize[2]'}
                rda.CSIMatrix_Size(3) = str2double(value);
            case {'PositionVector[0]'}
                rda.PositionVector(1) = str2double(value);
            case {'PositionVector[1]'}
                rda.PositionVector(2) = str2double(value);
            case {'PositionVector[2]'}
                rda.PositionVector(3) = str2double(value);
            case {'RowVector[0]'}
                rda.RowVector(1) = str2double(value);
            case {'RowVector[1]'}
                rda.RowVector(2) = str2double(value);
            case {'RowVector[2]'}
                rda.RowVector(3) = str2double(value);
            case {'ColumnVector[0]'}
                rda.ColumnVector(1) = str2double(value);
            case {'ColumnVector[1]'}
                rda.ColumnVector(2) = str2double(value);
            case {'ColumnVector[2]'}
                rda.ColumnVector(3) = str2double(value);
            otherwise
                % We don't know what this variable is.  Report this just to keep things clear
                %disp(['Unrecognised variable ' , variable ]);
        end

    else
        % Don't bother storing this bit of the output

    end

end

%%RE 110726 Take the used bits of the header info
if isfield(rda,'VOIRotationInPlane')
    MRS_struct.p.VoI_InPlaneRot(ii) = rda.VOIRotationInPlane;
else
    MRS_struct.p.VoI_InPlaneRot(ii) = 0;
end
MRS_struct.p.voxdim(ii,1) = rda.FoVHeight;
MRS_struct.p.voxdim(ii,2) = rda.FoVWidth;
MRS_struct.p.voxdim(ii,3) = rda.SliceThickness;
MRS_struct.p.voxoff(ii,1) = rda.PositionVector(1);
MRS_struct.p.voxoff(ii,2) = rda.PositionVector(2);
MRS_struct.p.voxoff(ii,3) = rda.PositionVector(3);
if isfield(rda,'MRFrequency')
    MRS_struct.p.LarmorFreq(ii) = rda.MRFrequency;
end
if isfield(rda,'VectorSize')
    MRS_struct.p.npoints(ii) = rda.VectorSize;
end
if isfield(rda,'DwellTime')
    MRS_struct.p.sw(ii) = 1/rda.DwellTime*1E6;
end
if isfield(rda,'TR')
    MRS_struct.p.TR(ii) = rda.TR;
end
if isfield(rda,'TE')
    MRS_struct.p.TE(ii) = rda.TE;
end

%
% So now we should have got to the point after the header text
%
% Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
%

% bytes_per_point = 16;
complex_data = fread(fid, rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.VectorSize * 2, 'double');

%fread(fid, 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
fclose(fid);

% Now convert this data into something meaningful

%Reshape so that we can get the real and imaginary separated
hmm = reshape(complex_data, 2, rda.VectorSize, rda.CSIMatrix_Size(1), rda.CSIMatrix_Size(2), rda.CSIMatrix_Size(3));

%Combine the real and imaginary into the complex matrix
hmm_complex = complex(hmm(1,:,:,:,:), hmm(2,:,:,:,:));

%RE 110726 This is the complex time domain data
MRS_struct.fids.offdata = hmm_complex;

%%%Now load in the ON data
rda_filename = on_filename; %This is generally file3
fid = fopen(rda_filename);

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';

tline = fgets(fid);

while isempty(strfind(tline, head_end_text))

    tline = fgets(fid);

    if isempty(strfind(tline, head_start_text)) + isempty(strfind(tline, head_end_text)) == 2

        % Store this data in the appropriate format
        occurence_of_colon = strfind(tline,':');
        variable = tline(1:occurence_of_colon-1);
        value    = tline(occurence_of_colon+1:length(tline));

        switch variable
            case {'PatientID', 'PatientName', 'StudyDescription', 'PatientBirthDate', 'StudyDate', 'StudyTime', 'PatientAge', 'SeriesDate', ...
                  'SeriesTime', 'SeriesDescription', 'ProtocolName', 'PatientPosition', 'ModelName', 'StationName', 'InstitutionName', ...
                  'DeviceSerialNumber', 'InstanceDate', 'InstanceTime', 'InstanceComments', 'SequenceName', 'SequenceDescription', 'Nucleus', ...
                  'TransmitCoil'}
                rda.(variable) = value;
            case {'PatientSex'}
                % Sex converter! (int to M,F,U)
                switch value
                    case 0
                        rda.sex = 'Unknown';
                    case 1
                        rda.sex = 'Male';
                    case 2
                        rda.sex = 'Female';
                end
            case {'SeriesNumber', 'InstanceNumber', 'AcquisitionNumber', 'NumOfPhaseEncodingSteps', 'NumberOfRows', 'NumberOfColumns', 'VectorSize'}
                %Integers
                rda.(variable) = str2double(value);
            case {'PatientWeight', 'TR', 'TE', 'TM', 'DwellTime', 'NumberOfAverages', 'MRFrequency', 'MagneticFieldStrength', 'FlipAngle', ...
                  'SliceThickness',  'FoVHeight', 'FoVWidth', 'PercentOfRectFoV', 'PixelSpacingRow', 'PixelSpacingCol'}
                %Floats
                rda.(variable) = str2double(value);
            case {'SoftwareVersion[0]'}
                rda.software_version = value;
            case {'CSIMatrixSize[0]'}
                rda.CSIMatrix_Size(1) = str2double(value);
            case {'CSIMatrixSize[1]'}
                rda.CSIMatrix_Size(2) = str2double(value);
            case {'CSIMatrixSize[2]'}
                rda.CSIMatrix_Size(3) = str2double(value);
            case {'PositionVector[0]'}
                rda.PositionVector(1) = str2double(value);
            case {'PositionVector[1]'}
                rda.PositionVector(2) = str2double(value);
            case {'PositionVector[2]'}
                rda.PositionVector(3) = str2double(value);
            case {'RowVector[0]'}
                rda.RowVector(1) = str2double(value);
            case {'RowVector[1]'}
                rda.RowVector(2) = str2double(value);
            case {'RowVector[2]'}
                rda.RowVector(3) = str2double(value);
            case {'ColumnVector[0]'}
                rda.ColumnVector(1) = str2double(value);
            case {'ColumnVector[1]'}
                rda.ColumnVector(2) = str2double(value);
            case {'ColumnVector[2]'}
                rda.ColumnVector(3) = str2double(value);
            otherwise
                % We don't know what this variable is.  Report this just to keep things clear
                %disp(['Unrecognised variable ' , variable ]);
        end

    else
        % Don't bother storing this bit of the output

    end

end

%%RE 110726 Take the used bits of the header info
MRS_struct.p.LarmorFreq(ii) = rda.MRFrequency;
MRS_struct.p.npoints(ii) = rda.VectorSize;
MRS_struct.p.Navg(ii) = rda.NumberOfAverages;
% So now we should have got to the point after the header text
%
% Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
%

% bytes_per_point = 16;
complex_data = fread(fid, rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.VectorSize * 2, 'double');

%fread(fid , 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
fclose(fid);

% Now convert this data into something meaningful

%Reshape so that we can get the real and imaginary separated
hmm = reshape(complex_data, 2, rda.VectorSize, rda.CSIMatrix_Size(1),  rda.CSIMatrix_Size(2),  rda.CSIMatrix_Size(3));

%Combine the real and imaginary into the complex matrix
hmm_complex = complex(hmm(1,:,:,:,:), hmm(2,:,:,:,:));

%RE 110726 This is the complex time domain data
MRS_struct.fids.ondata = hmm_complex;
MRS_struct.fids.data = [MRS_struct.fids.ondata; MRS_struct.fids.offdata].';

if nargin == 4

    %%%Now load in the Water data
    rda_filename=water_filename; %This is generally file3
    fid = fopen(rda_filename);

    head_start_text = '>>> Begin of header <<<';
    head_end_text   = '>>> End of header <<<';

    tline = fgets(fid);

    while (isempty(strfind(tline, head_end_text))) %#ok<*STREMP>

        tline = fgets(fid);

        if isempty(strfind(tline, head_start_text)) + isempty(strfind(tline, head_end_text)) == 2
            % Store this data in the appropriate format
            occurence_of_colon = strfind(':',tline);
            variable = tline(1:occurence_of_colon-1);
            value    = tline(occurence_of_colon+1:length(tline));

            switch variable
                case {'PatientID', 'PatientName', 'StudyDescription', 'PatientBirthDate', 'StudyDate', 'StudyTime', 'PatientAge', 'SeriesDate', ...
                      'SeriesTime', 'SeriesDescription', 'ProtocolName', 'PatientPosition', 'ModelName', 'StationName', 'InstitutionName', ...
                      'DeviceSerialNumber', 'InstanceDate', 'InstanceTime', 'InstanceComments', 'SequenceName', 'SequenceDescription', 'Nucleus', ...
                      'TransmitCoil'}
                    rda.(variable) = str2double(value);
                case {'PatientSex'}
                    % Sex converter! (int to M,F,U)
                    switch value
                        case 0
                            rda.sex = 'Unknown';
                        case 1
                            rda.sex = 'Male';
                        case 2
                            rda.sex = 'Female';
                    end
                case {'SeriesNumber', 'InstanceNumber', 'AcquisitionNumber', 'NumOfPhaseEncodingSteps', 'NumberOfRows', 'NumberOfColumns', 'VectorSize'}
                    %Integers
                    rda.(variable) = str2double(value);
                case {'PatientWeight', 'TR', 'TE', 'TM', 'DwellTime', 'NumberOfAverages', 'MRFrequency', 'MagneticFieldStrength', 'FlipAngle', ...
                      'SliceThickness', 'FoVHeight', 'FoVWidth', 'PercentOfRectFoV', 'PixelSpacingRow', 'PixelSpacingCol'}
                    %Floats
                    rda.(variable) = str2double(value);
                case {'SoftwareVersion[0]'}
                    rda.software_version = value;
                case {'CSIMatrixSize[0]'}
                    rda.CSIMatrix_Size(1) = str2double(value);
                case {'CSIMatrixSize[1]'}
                    rda.CSIMatrix_Size(2) = str2double(value);
                case {'CSIMatrixSize[2]'}
                    rda.CSIMatrix_Size(3) = str2double(value);
                case {'PositionVector[0]'}
                    rda.PositionVector(1) = str2double(value);
                case {'PositionVector[1]'}
                    rda.PositionVector(2) = str2double(value);
                case {'PositionVector[2]'}
                    rda.PositionVector(3) = str2double(value);
                case {'RowVector[0]'}
                    rda.RowVector(1) = str2double(value);
                case {'RowVector[1]'}
                    rda.RowVector(2) = str2double(value);
                case {'RowVector[2]'}
                    rda.RowVector(3) = str2double(value);
                case {'ColumnVector[0]'}
                    rda.ColumnVector(1) = str2double(value);
                case {'ColumnVector[1]'}
                    rda.ColumnVector(2) = str2double(value);
                case {'ColumnVector[2]'}
                    rda.ColumnVector(3) = str2double(value);
                otherwise
                    % We don't know what this variable is.  Report this just to keep things clear
                    %disp(['Unrecognised variable ' , variable ]);
            end

        else
            % Don't bother storing this bit of the output
        end

    end

    %%RE 110726 Take the used bits of the header info
    MRS_struct.p.LarmorFreq(ii) = rda.MRFrequency;
    MRS_struct.p.npoints(ii) = rda.VectorSize;
    %
    % So now we should have got to the point after the header text
    %
    % Siemens documentation suggests that the data should be in a double complex format (8bytes for real, and 8 for imaginary?)
    %

    %bytes_per_point = 16;
    complex_data = fread(fid, rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.CSIMatrix_Size(1) * rda.VectorSize * 2, 'double');

    %fread(fid, 1, 'double');  %This was a check to confirm that we had read all the data (it passed!)
    fclose(fid);

    % Now convert this data into something meaningful

    %Reshape so that we can get the real and imaginary separated
    hmm = reshape(complex_data, 2, rda.VectorSize, rda.CSIMatrix_Size(1), rda.CSIMatrix_Size(2), rda.CSIMatrix_Size(3));

    %Combine the real and imaginary into the complex matrix
    hmm_complex = complex(hmm(1,:,:,:,:), hmm(2,:,:,:,:));

    %RE 110726 This is the complex time domain data
    MRS_struct.fids.data_water = hmm_complex.';

end

end



