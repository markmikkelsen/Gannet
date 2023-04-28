function MRS_struct = SiemensRead(MRS_struct, off_filename, on_filename, water_filename)

ii = MRS_struct.ii;

% Load water-suppressed data
[off_data, hdr] = read_rda_data(off_filename);
on_data         = read_rda_data(on_filename);

MRS_struct.fids.data = [on_data; off_data].';

% Header info
if isfield(hdr, 'VOIRotationInPlane')
    MRS_struct.p.VoI_InPlaneRot(ii) = hdr.VOIRotationInPlane;
else
    MRS_struct.p.VoI_InPlaneRot(ii) = 0;
end

MRS_struct.p.LarmorFreq(ii) = hdr.MRFrequency;
MRS_struct.p.npoints(ii)    = hdr.VectorSize;
MRS_struct.p.sw(ii)         = 1/hdr.DwellTime * 1e6;
MRS_struct.p.TR(ii)         = hdr.TR;
MRS_struct.p.TE(ii)         = hdr.TE;
MRS_struct.p.Navg(ii)       = hdr.NumberOfAverages;
MRS_struct.p.voxdim(ii,1)   = hdr.FoVHeight;
MRS_struct.p.voxdim(ii,2)   = hdr.FoVWidth;
MRS_struct.p.voxdim(ii,3)   = hdr.SliceThickness;
MRS_struct.p.voxoff(ii,:)   = hdr.PositionVector;

% Load water data
if nargin == 4

    [water_data, hdr] = read_rda_data(water_filename);
    
    MRS_struct.fids.data_water = water_data.';

    if isfield(hdr, 'TR')
        MRS_struct.p.TR_water(ii) = hdr.TR;
    end

    if isfield(hdr, 'TE')
        MRS_struct.p.TE_water(ii) = hdr.TE;
    end

end

end


function [data, hdr] = read_rda_data(fname)

fid = fopen(fname);

head_start_text = '>>> Begin of header <<<';
head_end_text   = '>>> End of header <<<';

tline = fgets(fid);

while isempty(strfind(tline, head_end_text)) %#ok<*STREMP>

    tline = fgets(fid);

    if isempty(strfind(tline, head_start_text)) && isempty(strfind(tline, head_end_text))

        % Store this data in the appropriate format
        occurence_of_colon = strfind(tline,':');
        variable = tline(1:occurence_of_colon-1);
        value    = tline(occurence_of_colon+1:length(tline));

        switch variable

            case {'PatientID', 'PatientName', 'StudyDescription', 'PatientBirthDate', 'StudyDate', 'StudyTime', 'PatientAge', 'SeriesDate', ...
                  'SeriesTime', 'SeriesDescription', 'ProtocolName', 'PatientPosition', 'ModelName', 'StationName', 'InstitutionName', ...
                  'DeviceSerialNumber', 'InstanceDate', 'InstanceTime', 'InstanceComments', 'SequenceName', 'SequenceDescription', 'Nucleus', ...
                  'TransmitCoil'}
                hdr.(variable) = value;

            case 'PatientSex'
                % Sex converter (int to M,F,U)
                switch value
                    case 0
                        hdr.sex = 'Unknown';
                    case 1
                        hdr.sex = 'Male';
                    case 2
                        hdr.sex = 'Female';
                end

            case {'SeriesNumber', 'InstanceNumber', 'AcquisitionNumber', 'NumOfPhaseEncodingSteps', 'NumberOfRows', 'NumberOfColumns', 'VectorSize'}
                %Integers
                hdr.(variable) = str2double(value);

            case {'PatientWeight', 'TR', 'TE', 'TM', 'DwellTime', 'NumberOfAverages', 'MRFrequency', 'MagneticFieldStrength', 'FlipAngle', ...
                  'SliceThickness', 'FoVHeight', 'FoVWidth', 'PercentOfRectFoV', 'PixelSpacingRow', 'PixelSpacingCol', 'VOIRotationInPlane'}
                %Floats
                hdr.(variable) = str2double(value);

            case 'SoftwareVersion[0]'
                hdr.software_version = value;

            case 'CSIMatrixSize[0]'
                hdr.CSIMatrix_Size(1) = str2double(value);

            case 'CSIMatrixSize[1]'
                hdr.CSIMatrix_Size(2) = str2double(value);

            case 'CSIMatrixSize[2]'
                hdr.CSIMatrix_Size(3) = str2double(value);

            case 'PositionVector[0]'
                hdr.PositionVector(1) = str2double(value);

            case 'PositionVector[1]'
                hdr.PositionVector(2) = str2double(value);

            case 'PositionVector[2]'
                hdr.PositionVector(3) = str2double(value);

            case 'RowVector[0]'
                hdr.RowVector(1) = str2double(value);

            case 'RowVector[1]'
                hdr.RowVector(2) = str2double(value);

            case 'RowVector[2]'
                hdr.RowVector(3) = str2double(value);

            case 'ColumnVector[0]'
                hdr.ColumnVector(1) = str2double(value);

            case 'ColumnVector[1]'
                hdr.ColumnVector(2) = str2double(value);

            case 'ColumnVector[2]'
                hdr.ColumnVector(3) = str2double(value);

            case 'TablePosSag'
                hdr.TablePosition(1) = str2double(value);

            case 'TablePosCor'
                hdr.TablePosition(2) = str2double(value);

            case 'TablePosTra'
                hdr.TablePosition(3) = str2double(value);

        end

    end

end

data = fread(fid, hdr.CSIMatrix_Size(1) * hdr.CSIMatrix_Size(1) * hdr.CSIMatrix_Size(1) * hdr.VectorSize * 2, 'double');

% Reshape so that we can get the real and imaginary separated
data = reshape(data, 2, hdr.VectorSize, hdr.CSIMatrix_Size(1), hdr.CSIMatrix_Size(2), hdr.CSIMatrix_Size(3));

% Combine the real and imaginary into a complex matrix
data = complex(data(1,:,:,:,:), data(2,:,:,:,:));

fclose(fid);

end



