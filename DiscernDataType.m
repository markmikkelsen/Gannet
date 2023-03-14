function MRS_struct = DiscernDataType(filename, MRS_struct)

[~,~,ext] = fileparts(filename);

switch lower(ext)
    case '.7'
        MRS_struct.p.vendor = 'GE';
    case '.dat'
        MRS_struct.p.vendor = 'Siemens_twix';
    case '.data'
        MRS_struct.p.vendor = 'Philips_data';
    case '.dcm'
        MRS_struct.p.vendor = 'DICOM';
    case {'.gz','.nii'}
        MRS_struct.p.vendor = 'NIfTI';
    case '.ima'
        MRS_struct.p.vendor = 'Siemens_dicom';
    case '.raw'
        MRS_struct.p.vendor = 'Philips_raw';
    case '.rda'
        MRS_struct.p.vendor = 'Siemens_rda';
    case '.sdat'
        MRS_struct.p.vendor = 'Philips';
    otherwise
        error('Unrecognized file type! Extension must be .7, .dat, .data, .dcm, .gz, .ima, .nii, .raw, .rda, or, .sdat.');
end

end



