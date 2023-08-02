function MRS_struct = CoRegStandAlone(metabfile, struc)
% MRS_struct = CoRegStandAlone(metabfile, struc)
%   Reads the relevant header information from MRS files; performs
%   co-registration to a 3D structural image in either NIfTI (*.nii) (for
%   GE, Philips, and Siemens data) or DICOM (*.dcm) (for GE data) format
%   and creates a binary voxel mask; performs tissue segmentation using
%   SPM12; and returns the tissue fractions of gray matter, white matter,
%   and cerebrospinal fluid in the MRS voxel.
%
%   Requires:
%       - SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
%
%   Input:
%       metabfile - cell array containing the path(s) to an MRS data file
%                   in one of the following formats: Philips SDAT/SPAR,
%                   Siemens TWIX (*.dat), Siemens RDA (*.rda), GE (.7), or
%                   DICOM (*.ima, *.dcm).
%       struc     - cell array containing the path(s) to either a NIfTI
%                   file (*.nii) (for Philips and Siemens data) or a folder
%                   containing DICOMs (.*dcm) (for GE data).
%
%       It is possible to batch-process data. In this case, the number of
%       elements of the metabfile cell array needs to match the number of
%       elements of the struc cell array, with matching order.
%
%   Example:
%       MRS_struct = CoRegStandAlone({'MRS1.dat', 'MRS2.dat'}, {'nii.dat', 'nii.dat'})
%
%       This example will co-register and segment two Siemens MRS voxels to
%       the same structural.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-09-19)
%       goeltzs1@jhmi.edu
%   Updates:
%       Dr. Mark Mikkelsen (Weill Cornell Medicine)
%
%   History:
%       2018-09-19: First version of the code
%       2019-10-24: Minor bug fix (line 44)
%       2020-07-29: Some minor cosmetic changes
%       2022-06-03: Fixed bug related to target metabolite
%       2023-03-13: Update CSV filenames; prevent overwriting if CSV file
%                   already exists in the output directory

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loadFile = which('GannetLoad');
fileID = fopen(loadFile, 'rt');
str = fread(fileID, Inf, '*uchar');
fclose(fileID);
str = char(str(:)');
expression = '(?<field>MRS_struct.version.Gannet = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
MRS_struct.version.Gannet = out.version;

expression = '(?<field>MRS_struct.version.load   = )''(?<version>.*?)''';
out = regexp(str, expression, 'names');
MRS_struct.version.load = out.version;

MRS_struct.version.coregstandalone = '230729';

MRS_struct.ii = 0;
if size(metabfile,2) == 1
    metabfile = metabfile';
end
MRS_struct.metabfile  = metabfile;
MRS_struct.p.HERMES   = 0;
MRS_struct.p.numScans = length(metabfile);

% Flags
MRS_struct.p.mat    = 0; % Save results in *.mat file? (0 = NO, 1 = YES (default))
MRS_struct.p.csv    = 1; % Save results in *.csv file? (0 = NO, 1 = YES (default))
MRS_struct.p.vox    = {'vox1'}; % Name of the voxel
MRS_struct.p.target = {'GABAGlx'}; % Name of the target metabolite
MRS_struct.p.hide   = 0; % Do not display output figures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2. Determine data parameters from header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discern input data format
MRS_struct = DiscernDataType(metabfile{1}, MRS_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3. Load data from files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(metabfile) % Loop over all files in the batch (from metabfile)

    MRS_struct.ii = ii;

    switch MRS_struct.p.vendor

        case 'DICOM'
            MRS_struct = DICOMRead(MRS_struct, metabfile{ii});

        case 'GE'
            MRS_struct = GERead(MRS_struct, metabfile{ii});

        case 'NIfTI'
            MRS_struct = NIfTIMRSRead(MRS_struct, metabfile{ii});

        case 'Philips'
            MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});

        case 'Philips_data'
            MRS_struct = PhilipsDataRead(MRS_struct, metabfile{ii});

        case 'Philips_raw'
            MRS_struct = PhilipsRawRead(MRS_struct, metabfile{ii}, 3, 0);

        case 'Siemens_dicom'
            MRS_struct = SiemensDICOMRead(MRS_struct, metabfile{ii});

        case 'Siemens_rda'
            MRS_struct = SiemensRead(MRS_struct, metabfile{ii}, metabfile{ii});

        case 'Siemens_twix'
            MRS_struct = SiemensTWIXRead(MRS_struct, metabfile{ii});

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   4. Call coregister function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct = CoReg(MRS_struct, struc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   5. Call segment function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct = Seg(MRS_struct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   6. Clean up, save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save MRS_struct as mat file
MRS_struct = rmfield(MRS_struct,'fids');
if MRS_struct.p.mat
    save(fullfile(pwd, 'MRS_struct_CoRegStandAlone.mat'), 'MRS_struct', '-v7.3');
end

% Export MRS_struct fields into csv file
if MRS_struct.p.csv
    csv_name = fullfile(pwd, 'CoRegStandAlone_output.csv');
    if exist(csv_name, 'file')
        run_count = 1;
        csv_name = fullfile(pwd, ['CoRegStandAlone_output' num2str(run_count) '.csv']);
        while 1
            if exist(csv_name, 'file')
                run_count = run_count + 1;
                csv_name  = fullfile(pwd, ['CoRegStandAlone_output' num2str(run_count) '.csv']);
            else
                break
            end
        end
    end
    fprintf('\nExporting tissue segmentation results to %s...\n', csv_name);

    if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
        filename = MRS_struct.metabfile(:,1:2:end)';
    else
        filename = MRS_struct.metabfile';
    end
    for ii = 1:length(filename)
        [~,b,c] = fileparts(filename{ii});
        out.filename(ii,1) = cellstr([b c]);
    end

    out.fGM  = MRS_struct.out.vox1.tissue.fGM(:);
    out.fWM  = MRS_struct.out.vox1.tissue.fWM(:);
    out.fCSF = MRS_struct.out.vox1.tissue.fCSF(:);

    round2 = @(x) round(x*1e3)/1e3;
    T = table(out.filename, round2(out.fGM), round2(out.fWM), round2(out.fCSF), ...
        'VariableNames', {'Filename', 'fGM', 'fWM', 'fCSF'});
    writetable(T, csv_name);
end

end



