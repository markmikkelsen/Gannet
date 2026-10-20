function MRS_struct = GannetPreInitialise(MRS_struct)

% Acquisition parameters (0 = NO; 1 = YES)
    MRS_struct.p.target = {'GABAGlx'}; % Edited metabolite(s) of interest; permitted options are:
                                       % If MEGA-PRESS:
                                       %   {'GABA'}, {'GABAGlx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       % If HERMES:
                                       %   {'GABAGlx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
                                       % If HERCULES:
                                       %   {'GABAGlx','GSH'}
                                       % If phantom data:
                                       %   and MEGA-PRESS: {'GABA'}, {'Glx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       %   and HERMES: {'GABA','GSH'}, {'Glx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
    MRS_struct.p.HERMES       = 0;  % Data were acquired using HERMES
    MRS_struct.p.HERCULES     = 0;  % Data were acquired using HERCULES (if 1, MRS_struct.p.HERMES must be set to 1 as well)
    MRS_struct.p.PRIAM        = 0;  % Data were acquired using PRIAM
    MRS_struct.p.ON_OFF_order = []; % The editing order applied at acquisition.
                                    % If empty (the default: []), Gannet will determine the editing order automatically.
                                    % Otherwise, input 'offfirst' or 'onfirst' for MEGA-edited data;
                                    % or, for HERMES/HERCULES data, input a four-letter combination, such as 'ABCD' or 'CBAD', etc.
                                    % (see: ﻿doi:10.1016/j.neuroimage.2016.07.056)
    MRS_struct.p.seqorig = 'JHU'; % Origin of Philips MEGA-PRESS or GE HERMES sequences;
                                  % options are 'JHU' or 'Philips' if Philips, or 'Lythgoe' if GE (for HERMES only)
    MRS_struct.p.phantom = 0;     % Data are from a phantom (assumes phantom was scanned at room temperature)
    
% Analysis parameters (0 = NO; 1 = YES)
    MRS_struct.p.LB            = 3; % Exponential line-broadening (in Hz); default = 3 Hz; for phantom data, ~1.5 Hz is recommended
    MRS_struct.p.water_ECC     = 1; % Perform eddy current correction on water data
    MRS_struct.p.metab_ECC     = 0; % Perform eddy current correction on metabolite data (requires a water reference)
    MRS_struct.p.water_removal = 1; % Remove residual water signal in the DIFF spectrum using HSVD
    MRS_struct.p.alignment     = 'RobustSpecReg'; % Alignment method; options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                                  % 'Cr', 'Cho', 'NAA', 'H2O', or 'none' (recommended for phantom data)
    MRS_struct.p.use_prealign_ref   = 0; % In some cases, using RobustSpecReg to align HERMES/HERCULES data can result in
                                         % worse alignment compared to the pre-aligned data; setting this parameter to 1 will
                                         % make RobustSpecReg use the averaged pre-aligned subspectra as references to align the
                                         % averaged post-aligned subspectra, which may improve the final alignment
    MRS_struct.p.vox                = {'vox1'}; % For naming voxels, e.g., {'DLPFC'}; if data were acquired using PRIAM this could be,
                                                % e.g., {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.join               = 0; % Join multiple files (this can be performed in batch to join files across multiple subjects)
    MRS_struct.p.fit_resid_water    = 0; % Fit the residual water signal in the OFF spectrum to calculate a water suppression factor
    MRS_struct.p.weighted_averaging = 1; % Average subspectra using weighted averaging; otherwise, use arithmetic averaging
    MRS_struct.p.normalize          = 0; % Normalize voxel masks created by GannetCoRegister and GannetSegment to MNI space
                                         % and, if more than dataset has been run in the pipeline, create a mean overlap voxel
                                         % (note that this is only run if GannetSegment is run)
    
% Other options (0 = NO; 1 = YES)
    MRS_struct.p.mat = 0; % Save MRS_struct as a .mat file
    MRS_struct.p.csv = 0; % Extract useful results from the output structure MRS_struct and export them to a .csv file (applies to
                          % GannetFit, GannetSegment, and GannetQuantify only)
    MRS_struct.p.show_fits    = 1; % Show all peak fits of multi-peak models
    MRS_struct.p.append       = 1; % Append PDF outputs into one PDF (separately for each module) (requires export_fig in the Gannet
                                   % directory to be added to the search path and Ghostscript to be installed)
    MRS_struct.p.hide         = 1; % Do not display popup output figures
    MRS_struct.p.plot_subspec = 1; % Plot and print averaged ON/OFF or A–D subspectra (useful for debugging/QA)
    MRS_struct.p.debug        = 0; % Display detailed error reports and model fitting failure diagnostics

end
