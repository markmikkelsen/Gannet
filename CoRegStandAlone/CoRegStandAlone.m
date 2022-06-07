function MRS_struct = CoRegStandAlone(metabfile, struc)
% CoRegStandAlone(metabfile, struc)
%   Reads the relevant header information from MRS files; performs
%   co-registration to a 3D structural image in either NIfTI (*.nii) (for
%   Philips and Siemens data) or DICOM (*.dcm) (for GE data) format and
%   creates a binary voxel mask; performs tissue segmentation using SPM12;
%   and returns the tissue fractions of gray matter, white matter, and
%   cerebrospinal fluid in the MRS voxel.
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
%
%   History:
%       2018-09-19: First version of the code.
%       2019-10-24: Minor bug fix (line 44).
%       2020-07-29: Some minor cosmetic changes.
%       2022-06-03: Fixed bug related to target metabolite

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1. Pre-initialise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRS_struct.version.Gannet = '3.2.1';
MRS_struct.version.load = '220607';
MRS_struct.ii = 0;
if size(metabfile,2) == 1
    metabfile = metabfile';
end
MRS_struct.metabfile = metabfile;
MRS_struct.p.HERMES = 0;

% Flags
MRS_struct.p.mat = 1; % Save results in *.mat file? (0 = NO, 1 = YES (default)).
MRS_struct.p.csv = 1; % Save results in *.csv file? (0 = NO, 1 = YES (default)).
MRS_struct.p.vox = {'vox1'}; % Name of the voxel
MRS_struct.p.target = {'GABAGlx'}; % Name of the target metabolite

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

        case 'GE'
            MRS_struct = GERead(MRS_struct, metabfile{ii});

        case 'Siemens_twix'
            MRS_struct = SiemensTwixRead(MRS_struct, metabfile{ii});

        case 'Siemens_dicom'
            MRS_struct = SiemensDICOMRead(MRS_struct, metabfile{ii});

        case 'dicom'
            MRS_struct = DICOMRead(MRS_struct, metabfile{ii});

        case 'Siemens_rda'
            MRS_struct = SiemensRead(MRS_struct, metabfile{ii}, metabfile{ii});

        case 'Philips'
            MRS_struct = PhilipsRead(MRS_struct, metabfile{ii});

        case 'Philips_data'
            MRS_struct = PhilipsRead_data(MRS_struct, metabfile{ii});

        case 'Philips_raw'
            MRS_struct = PhilipsRawLoad(MRS_struct, metabfile{ii}, 3, 0);

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
    csv_name = fullfile(pwd, 'MRS_struct.csv');
    if exist(csv_name, 'file')
        fprintf('\nUpdating results in %s\n\n', 'MRS_struct.csv...');
    else
        fprintf('\nExporting results to %s\n\n', 'MRS_struct.csv...');
    end

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



