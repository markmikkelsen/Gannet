function MRS_struct = CalcIU(MRS_struct, vox, metab, ii)
% Function for quantifying concentration in institutional units
% Convert metabolite and water areas to institutional units
% (pseudo-concentration in mmol/L)

TR = MRS_struct.p.TR(ii)/1e3;
TE = MRS_struct.p.TE(ii)/1e3;
if isfield(MRS_struct.p,'TR_water')
    TR_water = MRS_struct.p.TR_water(ii)/1e3;
else
    TR_water = TR;
end
if isfield(MRS_struct.p,'TE_water')
    TE_water = MRS_struct.p.TE_water(ii)/1e3;
else
    TE_water = TE;
end
PureWaterConc   = 55.51*1e3; % mol/kg
WaterVisibility = 0.65; % this is approx the value from Ernst, Kreis, Ross (1993, JMR)
T1_Water        = 1.100; % average of WM and GM, Wansapura et al. 1999 (JMRI)
T2_Water        = 0.095; % average of WM and GM, Wansapura et al. 1999 (JMRI)
N_H_Water       = 2;

switch metab
    case 'GABA'
        EditingEfficiency = 0.5; % For TE = 68 ms
        T1_Metab  = 1.31;  % Puts et al. 2013 (JMRI)
        T2_Metab  = 0.088; % Edden et al. 2012 (JMRI)
        N_H_Metab = 2;
        MM        = 0.45; % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
                          % This fraction is platform- and implementation-dependent, based on length and
                          % shape of editing pulses and ifis Henry method
        
    case 'Glx'
        EditingEfficiency = 0.4; % determined by FID-A simulations (for TE = 68 ms)
        T1_Metab  = 1.23; % Posse et al. 2007 (MRM)
        T2_Metab  = 0.18; % Ganji et al. 2012 (NMR Biomed)
        N_H_Metab = 1;
        MM        = 1;
        
    case 'GSH'
        EditingEfficiency = 0.74; % Sanaei Nezhad et al., 2017, MRM, doi:10.1002/mrm.26532
        T1_Metab  = 0.397; % Choi et al., 2013, NMR Biomed., doi:10.1002/nbm.2815
        T2_Metab  = 0.088; % Choi et al., 2025, NMR Biomed., doi:10.1002/nbm.5313 
        N_H_Metab = 2;
        MM        = 1;
        
    case 'Lac'
        EditingEfficiency = 0.94; % determined by FID-A simulations (for TE = 140 ms)
        T1_Metab  = 1.50; % Wijnen et al. 2015 (NMR Biomed)
        T2_Metab  = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
        N_H_Metab = 3;
        MM        = 1;
        
    case 'EtOH'
        EditingEfficiency = 0.5; % assuming same as GABA for now
        T1_Metab  = 1.31;  % assuming same as GABA
        T2_Metab  = 0.088; % assuming same as GABA
        N_H_Metab = 3;
        MM        = 1;

    case 'Cr' % 3 ppm moiety
        EditingEfficiency = 1; % not edited, so 1
        T1_Metab  = (1.46 + 1.24)/2; % Mlynárik et al. 2001 (NMR in Biomed)
        T2_Metab  = (166 + 144 + 148)/3/1e3; % Wyss et al. 2018 (MRM)
        N_H_Metab = 3;
        MM        = 1;

    case 'Cho' % 3.2 ppm moiety
        EditingEfficiency = 1; % not edited, so 1
        T1_Metab  = (1.30 + 1.08)/2; % Mlynárik et al. 2001 (NMR in Biomed)
        T2_Metab  = (218 + 222 + 274)/3/1e3; % Wyss et al. 2018 (MRM)
        N_H_Metab = 9;
        MM        = 1;

    case 'NAA' % 2 ppm moiety
        EditingEfficiency = 1; % not edited, so 1
        T1_Metab  = (1.47 + 1.35)/2; % Mlynárik et al. 2001 (NMR in Biomed)
        T2_Metab  = (343 + 263 + 253)/3/1e3; % Wyss et al. 2018 (MRM)
        N_H_Metab = 3;
        MM        = 1;

    case 'Glu' % 2.34 ppm moiety
        EditingEfficiency = 0.4; % Saleh et al. 2024 (MRM)
        T1_Metab  = 1.23; % Posse et al. 2007 (MRM)
        T2_Metab  = 0.18; % Ganji et al. 2012 (NMR Biomed)
        N_H_Metab = 2;
        MM        = 1;
end

T1_Factor = (1 - exp(-TR_water./T1_Water)) ./ (1 - exp(-TR./T1_Metab));
T2_Factor = exp(-TE_water./T2_Water) ./ exp(-TE./T2_Metab);

if strcmp(MRS_struct.p.vendor, 'Siemens_rda')
    % Factor of 2 is appropriate for averaged Siemens data (read in separately as ON and OFF)
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii)) ...
        .* PureWaterConc .* WaterVisibility .* T1_Factor .* T2_Factor .* (N_H_Water ./ N_H_Metab) ...
        .* MM ./ 2 ./ EditingEfficiency;
else
    MRS_struct.out.(vox).(metab).ConcIU(ii) = (MRS_struct.out.(vox).(metab).Area(ii) ./ MRS_struct.out.(vox).water.Area(ii)) ...
        .* PureWaterConc .* WaterVisibility .* T1_Factor .* T2_Factor .* (N_H_Water ./ N_H_Metab) ...
        .* MM ./ EditingEfficiency;
end



