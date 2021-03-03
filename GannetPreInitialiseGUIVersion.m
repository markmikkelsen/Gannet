function MRS_struct = GannetPreInitialiseGUIVersion(configPath, MRS_struct)

% Load the given configuration file
load(configPath, 'selectedMetabolites', 'orderOfEditingPulses', 'philipsPatch',...
    'alignmentSelection', 'lineBroadening', 'eddyCorrectionOnMetabolite',...
    'eddyCorrectionOnWater','removeResidualWaterSignal','fitResidualWaterSignal',...
    'voxel1Name','voxel2Name','saveProcessDifferenceSpectrum','extractUsefulDataFromMRSStructure',...
    'saveMrsStructAsMat','hermesFlag', 'herculesFlag', 'priamFlag', 'phantomFlag', 'preAlignmentFlag',...
    'weightedAveragingFlag', 'silentFlag');

% Acquisition parameters
    MRS_struct.p.target = selectedMetabolites; % edited metabolite(s) of interest; allowable options are:
                                       % if MEGA-PRESS:
                                       %   {'GABAGlx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       % if HERMES:
                                       %   {'GABAGlx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
                                       % if HERCULES:
                                       %   {'GABAGlx','GSH'}
                                       % if phantom data:
                                       %   and MEGA-PRESS: {'GABA'}, {'Glx'}, {'GSH'}, {'Lac'}, or {'EtOH'}
                                       %   and HERMES: {'GABA','GSH'}, {'Glx','GSH'}, {'Lac','GSH'}, or {'EtOH','GABA','GSH'}
    MRS_struct.p.ON_OFF_order = orderOfEditingPulses; % order of editing pulses; options are 'onfirst' or 'offfirst'
    MRS_struct.p.seqorig      = philipsPatch; % origin of Philips MEGA-PRESS or GE HERMES sequences;
                                     % options are 'JHU' or 'Philips' if Philips, or 'Lythgoe' or 'Noeske' if GE (HERMES only)
    
% Analysis parameters
    MRS_struct.p.LB            = lineBroadening; % exponential line-broadening (in Hz)
    MRS_struct.p.water_ECC     = eddyCorrectionOnWater; % 1 = YES, perform eddy current correction on water data
    MRS_struct.p.metab_ECC     = eddyCorrectionOnMetabolite; % 1 = YES, perform eddy current correction on metabolite data
    MRS_struct.p.water_removal = removeResidualWaterSignal; % 1 = YES, remove residual water signal using HSVD
    MRS_struct.p.alignment     = alignmentSelection; % alignment method; options are 'RobustSpecReg' (recommended), 'SpecReg', 'SpecRegHERMES',
                                                  % 'Cr', 'Cho', 'NAA', 'H2O', 'CrOFF', or 'none' (recommended for phantom data)
    MRS_struct.p.use_prealign_ref = preAlignmentFlag; % 1 = YES; in some cases, using RobustSpecReg to align HERMES/HERCULES data can result in
                                       % worse alignment compared to the pre-aligned data; setting this parameter to 1 will
                                       % make RobustSpecReg use the averaged pre-aligned subspectra as references to align the
                                       % averaged post-aligned subspectra, which may improve the final alignment
    MRS_struct.p.vox = {voxel1Name, voxel2Name}; % for naming voxels in PRIAM data, e.g. {'anterior','posterior'}, {'right','left'}, etc.
    MRS_struct.p.fit_resid_water = fitResidualWaterSignal; % 1 = YES, fit the residual water signal in the DIFF spectrum to calculate water suppression factor
    MRS_struct.p.weighted_averaging = weightedAveragingFlag; % 1 = YES, average subspectra using weighted averaging
    
% Flags
    MRS_struct.p.HERMES   = hermesFlag; % 1 = YES, 0 = NO
    MRS_struct.p.HERCULES = herculesFlag; % 1 = YES, 0 = NO (if 1, MRS_struct.p.HERMES *must* be set to 1 as well)
    MRS_struct.p.PRIAM    = priamFlag; % 1 = YES, 0 = NO
    MRS_struct.p.phantom  = phantomFlag; % 1 = YES (assumes phantom was scanned at room temperature), 0 = NO (for in vivo data)
    MRS_struct.p.mat      = saveMrsStructAsMat; % 1 = YES, save MRS_struct as .mat file
    MRS_struct.p.sdat     = saveProcessDifferenceSpectrum; % 1 = YES, save processed DIFF spectrum as .sdat file (only for Philips SDAT MEGA-PRESS data)
    MRS_struct.p.csv      = extractUsefulDataFromMRSStructure; % 1 = YES, extract useful data from MRS_struct and export to .csv file (applies to GannetFit,
                               % GannetSegment and GannetQuantify)
    MRS_struct.p.silent   = silentFlag; % 1 = YES, do not dynamically display output figures
    
end



