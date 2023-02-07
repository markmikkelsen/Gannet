function MRS_struct = GannetPreInitialise(MRS_struct)

% Acquisition parameters
    MRS_struct.p.target = {'Lac'}; % Edited metabolite(s) of interest; permitted options are:
                                       % If MEGA-PRESS:
                                       %   {'GABA'}, {'GABAGlx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       % If HERMES:
                                       %   {'GABAGlx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
                                       % If HERCULES:
                                       %   {'GABAGlx','GSH'}
                                       % If phantom data:
                                       %   and MEGA-PRESS: {'GABA'}, {'Glx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       %   and HERMES: {'GABA','GSH'}, {'Glx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
    MRS_struct.p.seqorig = 'JHU'; % Origin of Philips MEGA-PRESS or GE HERMES sequences;
                                  % options are 'JHU' or 'Philips' if Philips, or 'Lythgoe' if GE (HERMES only)
    
% Analysis parameters
    MRS_struct.p.LB            = 3; % Exponential line-broadening (in Hz)
    MRS_struct.p.water_ECC     = 1; % 1 = YES, perform eddy current correction on water data
    MRS_struct.p.metab_ECC     = 0; % 1 = YES, perform eddy current correction on metabolite data (requires a water reference)
    MRS_struct.p.water_removal = 1; % 1 = YES, remove residual water signal in DIFF spectrum using HSVD
    MRS_struct.p.alignment     = 'RobustSpecReg'; % Alignment method; options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                                  % 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF', or 'none' (recommended for phantom data)
    MRS_struct.p.use_prealign_ref = 0; % 1 = YES; in some cases, using RobustSpecReg to align HERMES/HERCULES data can result in
                                       % worse alignment compared to the pre-aligned data; setting this parameter to 1 will
                                       % make RobustSpecReg use the averaged pre-aligned subspectra as references to align the
                                       % averaged post-aligned subspectra, which may improve the final alignment
    MRS_struct.p.vox                = {'vox1'}; % For naming voxels in PRIAM data, e.g. {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.fit_resid_water    = 0; % 1 = YES, fit the residual water signal in the OFF spectrum to calculate a water suppression factor
    MRS_struct.p.weighted_averaging = 1; % 1 = YES, average subspectra using weighted averaging
    
% Flags(0 = NO; 1 = YES)
    MRS_struct.p.HERMES   = 0; % Data acquired using HERMES
    MRS_struct.p.HERCULES = 0; % Data acquired using HERCULES; if 1, MRS_struct.p.HERMES must be set to 1 as well
    MRS_struct.p.PRIAM    = 0; % Data acquired using PRIAM
    MRS_struct.p.phantom  = 0; % Data are from a phantom (assumes phantom was scanned at room temperature)
    MRS_struct.p.join     = 0; % Join multiple files (this can be batched across subjects)
    MRS_struct.p.mat      = 0; % Save MRS_struct as a .mat file
    MRS_struct.p.csv      = 0; % Extract useful data from MRS_struct and export them to a .csv file (applies to GannetFit,
                               % GannetSegment and GannetQuantify)
    MRS_struct.p.append   = 1; % Append PDF outputs into one PDF (separately for each module) (requires export_fig in the Gannet
                               % directory to be added to the search path and Ghostscript to be installed)
    MRS_struct.p.hide     = 1; % Do not display output figures
    
end



