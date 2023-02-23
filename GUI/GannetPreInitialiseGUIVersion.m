function MRS_struct = GannetPreInitialiseGUIVersion(configFilePath, MRS_struct)

try
    load(configFilePath, 'configFile');
    
    % Acquisition parameters
    MRS_struct.p.target = configFile.metabolitesOfInterest; % edited metabolite(s) of interest; allowable options are:
                                                            % if MEGA-PRESS:
                                                            %   {'GABAGlx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                                            % if HERMES:
                                                            %   {'GABAGlx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
                                                            % if HERCULES:
                                                            %   {'GABAGlx','GSH'}
                                                            % if phantom data:
                                                            %   and MEGA-PRESS: {'GABA'}, {'Glx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                                            %   and HERMES: {'GABA','GSH'}, {'Glx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
    MRS_struct.p.seqorig = configFile.originOfSequences; % origin of Philips MEGA-PRESS or GE HERMES sequences;
                                                         % options are 'JHU' or 'Philips' if Philips, or 'Lythgoe' if GE (HERMES only)
    
    % Analysis parameters
    MRS_struct.p.LB            = configFile.lineBroadening; % exponential line-broadening (in Hz)
    MRS_struct.p.water_ECC     = configFile.ECCWater; % 1 = YES, perform eddy current correction on water data
    MRS_struct.p.metab_ECC     = configFile.ECCMetab; % 1 = YES, perform eddy current correction on metabolite data
    MRS_struct.p.water_removal = configFile.removeResidualWater; % 1 = YES, remove residual water signal using HSVD
    MRS_struct.p.alignment     = configFile.alignmentMethod; % alignment method; options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                                             % 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF', or 'none' (recommended for phantom data)
    MRS_struct.p.use_prealign_ref = configFile.usePreAlignedRef; % 1 = YES; in some cases, using RobustSpecReg to align HERMES/HERCULES data can result in
                                                             % worse alignment compared to the pre-aligned data; setting this parameter to 1 will
                                                             % make RobustSpecReg use the averaged pre-aligned subspectra as references to align the
                                                             % averaged post-aligned subspectra, which may improve the final alignment
    MRS_struct.p.vox                = configFile.voxelNames; % for naming voxels in PRIAM data, e.g. {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.fit_resid_water    = configFile.fitResidualWater; % 1 = YES, fit the residual water signal in the DIFF spectrum to calculate water suppression factor
    MRS_struct.p.weighted_averaging = configFile.weightedAveraging; % 1 = YES, average subspectra using weighted averaging
    
    % Flags
    MRS_struct.p.HERMES   = configFile.isHERMES;    % 1 = YES, 0 = NO
    MRS_struct.p.HERCULES = configFile.isHERCULES;  % 1 = YES, 0 = NO (if 1, MRS_struct.p.HERMES *must* be set to 1 as well)
    MRS_struct.p.PRIAM    = configFile.isPRIAM;     % 1 = YES, 0 = NO
    MRS_struct.p.phantom  = configFile.isPhantom;   % 1 = YES (assumes phantom was scanned at room temperature), 0 = NO (for in vivo data)
    MRS_struct.p.join     = configFile.join;        % Join multiple files (this can be batched across subjects)
    MRS_struct.p.mat      = configFile.saveAsMAT;   % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.csv      = configFile.exportToCSV; % 1 = YES, extract useful data from MRS_struct and export to .csv file (applies to GannetFit, GannetSegment and GannetQuantify)
    MRS_struct.p.append   = configFile.append;      % Append PDF outputs into one PDF (separately for each module) (requires export_fig in the Gannet
                                                    % folder to be added to the search path and GhostScript to be installed)
    MRS_struct.p.hide     = configFile.hide;        % Do not display output figures
    
    MRS_struct.p.readWriteAccess    = configFile.readWriteAccess; % 0 = Missing read/write access to input data
    MRS_struct.p.errorReadingConfig = false; % Set flag to false if config was read properly; stop GannetLoad if there was an issue
catch
    MRS_struct.p.errorReadingConfig = true;
    warning(strcat("Error reading configuration file: ", configFilePath, ". Make sure all parameters are correct. Stopping GannetLoad."));
end
end



