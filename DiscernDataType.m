function MRS_struct = DiscernDataType(filename, MRS_struct)

[~,~,ext] = fileparts(filename);

switch lower(ext)
    case '.7'
        MRS_struct.p.vendor = 'GE';
    case '.sdat'
        MRS_struct.p.vendor = 'Philips';
        if all(isstrprop(ext(end-3:end), 'upper'))
            MRS_struct.p.spar_string = 'SPAR';
        else
            MRS_struct.p.spar_string = 'spar';
        end
    case '.data'
        MRS_struct.p.vendor = 'Philips_data';
    case '.raw'
        MRS_struct.p.vendor = 'Philips_raw';
    case '.rda'
        MRS_struct.p.vendor = 'Siemens_rda';
    case '.dat'
        MRS_struct.p.vendor = 'Siemens_twix';
    case '.ima'
        MRS_struct.p.vendor = 'Siemens_dicom';
    case '.dcm'
        MRS_struct.p.vendor = 'dicom';
    otherwise
        error('Unrecognized file type! Extension should be .7, .sdat, .data, .raw, .rda, .dat, .ima, or .dcm.');
end

end



