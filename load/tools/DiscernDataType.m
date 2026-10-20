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

if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
    w = warning('query','backtrace');
    if strcmp(w.state,'on')
        warning('off','backtrace');
    end
    fprintf('\n');
    warning(['The Siemens .rda format is NOT recommended for use in Gannet. ' ...
             'If possible, please re-export your data in the TWIX (.dat) format.']);
    if strcmp(w.state,'on')
        warning('on','backtrace');
    end
end

end



