function DicomHeader = read_dcm_header(fid)
% DicomHeader = read_dcm_header(fid)
%   Reads header information from a DICOM file.
%
%   Example:
%       dcmHeader = read_dcm_header('file0001.dcm')
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-04-24)
%       goeltzs1@jhmi.edu
%
%   Credits:
%
%   Version history:
%   0.9:  First version (2018-04-24)
%   0.91: Several sequence-specific loading fixes (2018-05-13)
%   0.92: Added support for sLASER sequence (2018-07-18)
%   0.93: Bug fix for invalid field names (thanks to Meredith Reid) (2022-04-27)
%   0.94: Added CMRR PRESS sequence type (2022-10-13)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HEADER INFO PARSING %%%
% Simply open the dicom file, the information should all be the same.
fid = fopen(fid);

% Start looking for a convenient parameter block. The line before will
% start with ### ASCCONV BEGIN and end with ### ASCCONV END.
% This is defined here.
head_start_text = '### ASCCONV BEGIN';
head_end_text   = '### ASCCONV END';
tline = fgets(fid); % get first line

% Keep looking until start of the parameter block is found.
while isempty(strfind(tline, head_start_text)) %#ok<*STREMP>
    tline = fgets(fid);
end

% Look for regular expression containing the 'equal' signs
while isempty(strfind(tline, head_end_text))
    [tokens, ~] = regexp(tline,'([\w\[\].]*)\s*=\s*([\w.-\"\\]*)','tokens','match');
    % When a matching string is found, parse the results into a struct
    if isscalar(tokens)
        fieldname = regexprep(tokens{1}{1}, '\[|\]|_',''); % delete invalid characters
        if isempty(strfind(tokens{1}{2},'"'))
            if strcmp(tokens{1}{2},'0x1')
                value = 0;
            elseif strcmp(tokens{1}{2},'1x1')
                value = 1;
            else
                value = str2double(tokens{1}{2});
            end
        else
            value = tokens{1}{2};
        end
        C = strsplit(fieldname,'.'); % check for nested variable names
        switch length(C)
            case 1
                try
                    dcmHeader.(C{1}) = value;
                catch
                    tline = fgets(fid);
                    continue;
                end
            case 2
                try
                    dcmHeader.(C{1}).(C{2}) = value;
                catch
                    tline = fgets(fid);
                    continue;
                end
            case 3
                dcmHeader.(C{1}).(C{2}).(C{3}) = value;
            case 4
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}) = value;
            case 5
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}) = value;
            case 6
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}).(C{6}) = value;
            case 7
                dcmHeader.(C{1}).(C{2}).(C{3}).(C{4}).(C{5}).(C{6}).(C{7}) = value;
        end
    end
    tline = fgets(fid);
end
fclose(fid);


% Determine sequence name and type
DicomHeader.sequenceFileName = dcmHeader.tSequenceFileName; % Full sequence name
% Determine the origin of the sequence
if strfind(DicomHeader.sequenceFileName,'svs_edit')
    DicomHeader.seqtype = 'MEGAPRESS';
    if strfind(DicomHeader.sequenceFileName(end-3:end),'univ') %#ok<STRIFCND>
        DicomHeader.seqorig = 'Universal'; % Universal sequence
    else
        DicomHeader.seqorig = 'WIP'; % Siemens WIP
    end
elseif strfind(DicomHeader.sequenceFileName,'JEdit')
    DicomHeader.seqtype = 'MEGAPRESS';
    DicomHeader.seqorig = 'Utah'; % Utah sequence
elseif strfind(DicomHeader.sequenceFileName,'jn_')
    DicomHeader.seqtype = 'MEGAPRESS';
    DicomHeader.seqorig = 'JN'; % Jamie Near's sequence
elseif strfind(DicomHeader.sequenceFileName,'eja_svs_mpress')
    DicomHeader.seqtype = 'MEGAPRESS';
    DicomHeader.seqorig = 'CMRR'; % Minnesota sequence
elseif strfind(DicomHeader.sequenceFileName,'eja_svs_press')
    DicomHeader.seqtype = 'PRESS';
    DicomHeader.seqorig = 'CMRR'; % Minnesota sequence
elseif strfind(DicomHeader.sequenceFileName,'svs_se')
    DicomHeader.seqtype = 'PRESS'; % PRESS
elseif strfind(DicomHeader.sequenceFileName,'svs_slaser')
    DicomHeader.seqtype = 'sLASER'; % sLASER
elseif strfind(DicomHeader.sequenceFileName,'st_vapor_643')
    DicomHeader.seqtype = 'STEAM'; % STEAM
elseif strfind(DicomHeader.sequenceFileName, '%SiemensSeq%\svs_st')
    DicomHeader.seqtype = 'STEAM'; % Siemens 7T Terra.X product STEAM sequence
else
    DicomHeader.seqorig = DicomHeader.sequenceFileName;
    error(['Unknown sequence: ' DicomHeader.seqorig '. Please consult the Gannet team for support.'])
end

% Read information
DicomHeader.TR = dcmHeader.alTR0 * 1e-3; % TR [ms]
DicomHeader.TE = dcmHeader.alTE0 * 1e-3; % TE [ms]
if isfield(dcmHeader, 'lAverages')
    DicomHeader.nAverages = dcmHeader.lAverages;
else
    % Minnesota sequence (CMRR, Eddy Auerbach) may store numbers of averages in a
    % different field. GO 112017. Spelling may vary as well...
    if isfield(dcmHeader, 'sWipMemBlock')
        DicomHeader.nAverages = dcmHeader.sWipMemBlock.alFree2;
    elseif isfield(dcmHeader, 'sWiPMemBlock')
        DicomHeader.nAverages = dcmHeader.sWiPMemBlock.alFree2;
    end
end
DicomHeader.removeOS   = dcmHeader.sSpecPara.ucRemoveOversampling; % Is the oversampling removed in the RDA files?
DicomHeader.vectorSize = dcmHeader.sSpecPara.lVectorSize; % Data points specified on exam card
% GO180424: If a parameter is set to zero (e.g. if no voxel rotation is
% performed), the respective field does not show up in the dicom file. This
% case needs to be intercepted. Setting to the minimum possible value.
if ~isfield(dcmHeader.sSpecPara.sVoI, 'dInPlaneRot')
    dcmHeader.sSpecPara.sVoI.dInPlaneRot = realmin('double');
end
VoI_Params = {'dCor','dSag','dTra'};
for pp = 1:length(VoI_Params)
    if ~isfield(dcmHeader.sSpecPara.sVoI.sNormal, VoI_Params{pp})
        dcmHeader.sSpecPara.sVoI.sNormal.(VoI_Params{pp}) = realmin('double');
    end
    if isfield(dcmHeader.sSpecPara.sVoI, 'sPosition')
        if ~isfield(dcmHeader.sSpecPara.sVoI.sPosition, VoI_Params{pp})
            dcmHeader.sSpecPara.sVoI.sPosition.(VoI_Params{pp}) = realmin('double');
        end
    end
end

DicomHeader.VoI_InPlaneRot = dcmHeader.sSpecPara.sVoI.dInPlaneRot; % Voxel rotation in plane
DicomHeader.VoI_RoFOV      = dcmHeader.sSpecPara.sVoI.dReadoutFOV; % Voxel size in readout direction [mm]
DicomHeader.VoI_PeFOV      = dcmHeader.sSpecPara.sVoI.dPhaseFOV; % Voxel size in phase encoding direction [mm]
DicomHeader.VoIThickness   = dcmHeader.sSpecPara.sVoI.dThickness; % Voxel size in slice selection direction [mm]
DicomHeader.NormCor        = dcmHeader.sSpecPara.sVoI.sNormal.dCor; % Coronal component of normal vector of voxel
DicomHeader.NormSag        = dcmHeader.sSpecPara.sVoI.sNormal.dSag; % Sagittal component of normal vector of voxel
DicomHeader.NormTra        = dcmHeader.sSpecPara.sVoI.sNormal.dTra; % Transversal component of normal vector of voxel
if isfield(dcmHeader.sSpecPara.sVoI, 'sPosition')
    DicomHeader.PosCor         = dcmHeader.sSpecPara.sVoI.sPosition.dCor; % Coronal coordinate of voxel [mm]
    DicomHeader.PosSag         = dcmHeader.sSpecPara.sVoI.sPosition.dSag; % Sagittal coordinate of voxel [mm]
    DicomHeader.PosTra         = dcmHeader.sSpecPara.sVoI.sPosition.dTra; % Transversal coordinate of voxel [mm]
end
% delta frequency (center of slice selection)
if isfield(dcmHeader.sSpecPara, 'dDeltaFrequency')
    DicomHeader.deltaFreq = dcmHeader.sSpecPara.dDeltaFrequency;
else
    DicomHeader.deltaFreq = 0;
end

if isfield(dcmHeader.sProtConsistencyInfo, 'flNominalB0')
    DicomHeader.B0 = dcmHeader.sProtConsistencyInfo.flNominalB0; % Nominal B0 [T]
end
DicomHeader.dwellTime = dcmHeader.sRXSPEC.alDwellTime0; % Dwell time [ns]
DicomHeader.tx_freq   = dcmHeader.sTXSPEC.asNucleusInfo0.lFrequency; % Transmitter frequency [Hz]

% These may only be extractable from a few sequences and MEGA-PRESS
% versions:
% Editing pulse parameters
if isfield(DicomHeader, 'seqorig')
    if strcmp(DicomHeader.seqorig, 'CMRR')
        if isfield(dcmHeader, 'sWipMemBlock')
            if isfield(dcmHeader.sWipMemBlock, 'adFree3')
                DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree3;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree2')
                DicomHeader.editRF.freq(2) = dcmHeader.sWipMemBlock.adFree2;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree6')
                DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree8;
            end
        elseif isfield(dcmHeader, 'sWiPMemBlock')
            if isfield(dcmHeader.sWiPMemBlock, 'adFree3')
                DicomHeader.editRF.freq(1) = dcmHeader.sWiPMemBlock.adFree3;
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree2')
                DicomHeader.editRF.freq(2) = dcmHeader.sWiPMemBlock.adFree2;
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree6')
                DicomHeader.editRF.bw = dcmHeader.sWiPMemBlock.adFree8;
            end
        end
    elseif strcmp(DicomHeader.seqorig, 'Utah')
        if isfield(dcmHeader, 'sWipMemBlock')
            if isfield(dcmHeader.sWipMemBlock, 'adFree8')
                DicomHeader.editRF.centerFreq = dcmHeader.sWipMemBlock.adFree8;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree3')
                DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree3;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree4')
                DicomHeader.editRF.freq(2) = dcmHeader.sWipMemBlock.adFree4;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree0')
                DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree0;
            end
        elseif isfield(dcmHeader, 'sWiPMemBlock')
            if isfield(dcmHeader.sWipMemBlock, 'adFree8')
                DicomHeader.editRF.centerFreq = dcmHeader.sWipMemBlock.adFree8;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree3')
                DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree3;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree4')
                DicomHeader.editRF.freq(2) = dcmHeader.sWipMemBlock.adFree4;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree0')
                DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree0;
            end
        end
    else
        if isfield(dcmHeader, 'sWipMemBlock')
            if isfield(dcmHeader.sWipMemBlock, 'adFree9')
                DicomHeader.editRF.centerFreq = dcmHeader.sWipMemBlock.adFree9;
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree7')
                DicomHeader.editRF.freq(1) = dcmHeader.sWipMemBlock.adFree7;
                DicomHeader.editRF.freq(2) = DicomHeader.editRF.centerFreq + (DicomHeader.editRF.centerFreq - DicomHeader.editRF.freq(1));
            end
            if isfield(dcmHeader.sWipMemBlock, 'adFree8')
                DicomHeader.editRF.bw = dcmHeader.sWipMemBlock.adFree8;
            end
        elseif isfield(dcmHeader, 'sWiPMemBlock')
            if isfield(dcmHeader.sWiPMemBlock, 'adFree9')
                DicomHeader.editRF.centerFreq = dcmHeader.sWiPMemBlock.adFree9;
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree7')
                DicomHeader.editRF.freq(1) = dcmHeader.sWiPMemBlock.adFree7;
                DicomHeader.editRF.freq(2) = DicomHeader.editRF.centerFreq + (DicomHeader.editRF.centerFreq - DicomHeader.editRF.freq(1));
            end
            if isfield(dcmHeader.sWiPMemBlock, 'adFree8')
                DicomHeader.editRF.bw = dcmHeader.sWiPMemBlock.adFree8;
            end
        end
    end
end



