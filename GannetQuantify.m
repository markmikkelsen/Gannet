function MRS_struct = GannetQuantify(MRS_struct)
% Function to derive absolute concentrations of metabolite signal estimates

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

MRS_struct.version.quantify = '230621';

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end

run_count = 0;

% Check if there are water files, otherwise exit
if ~strcmp(MRS_struct.p.reference, 'H2O')
    fprintf('\n');
    error('No water reference files found in input structure ''%s''. GannetQuantify.m requires water references to run. Exiting...', inputname(1));
end

% ******
% RAEE (190107): Major change to water concentration calc to bring into
% line with Gasparovic et al. 2006. Tagged '% Gasparovic et al. method
% (RAEE)'. Fractions implemented as molar fractions, not volume fractions.
% ******

% Constants
% From Wansapura et al. 1999 (JMRI)
%        T1          T2
% WM   832 +/- 10  79.2 +/- 0.6
% GM  1331 +/- 13  110 +/- 2
%
% From Lu et al. 2005 (JMRI)
% CSF T1 = 3817 +/- 424msec - but state may underestimated and that 4300ms
% is likely more accurate - but the reference is to an ISMRM 2001 abstract
% MacKay (last author) 2006 ISMRM abstract has T1 CSF = 3300 ms
% CSF T2 = 503.0 +/- 64.3 Piechnik MRM 2009; 61: 579
% However, other values from Stanisz et al:
% CPMG for T2, IR for T1
% T2GM = 99 +/ 7, lit: 71 +/- 10 (27)
% T1GM = 1820 +/- 114, lit 1470 +/- 50 (29)
% T2WM = 69 +/-3 lit 56 +/- 4 (27)
% T1WM = 1084 +/- 45 lit 1110 +/- 45 (29)

T1w_WM    = 0.832;
T2w_WM    = 0.0792;
T1w_GM    = 1.331;
T2w_GM    = 0.110;
T1w_CSF   = 3.817;
T2w_CSF   = 0.503;
N_H_Water = 2;

% Determine concentration of water in GM, WM and CSF
% Gasparovic et al. 2006 (MRM) uses relative densities,
% ref to Ernst et al. 1993 (JMR)
% fGM = 0.78
% fWM = 0.65
% fCSF = 0.97
% such that
% concw_GM  = 0.78 * 55.51 mol/kg = 43.30
% concw_WM  = 0.65 * 55.51 mol/kg = 36.08
% concw_CSF = 0.97 * 55.51 mol/kg = 53.84

concW_GM    = 43.30*1e3;
concW_WM    = 36.08*1e3;
concW_CSF   = 53.84*1e3;
molal_concW = 55.51*1e3; % Gasparovic et al. method (RAEE)

% Loop over voxels if PRIAM
for kk = 1:length(vox)

    meanfGM = mean(MRS_struct.out.(vox{kk}).tissue.fGM); % average GM fractions across subjects
    meanfWM = mean(MRS_struct.out.(vox{kk}).tissue.fWM); % average WM fractions across subjects

    for ii = 1:MRS_struct.p.numScans

        if kk == 1 && ii == 1
            fprintf('\nQuantifying metabolites...\n');
        end

        target = [MRS_struct.p.target, {'Cr'}, {'Cho'}, {'NAA'}]; % Add Cr, Cho, and NAA
        tmp    = strcmp(target,'GABAGlx');
        if any(tmp)
            target = {'GABA','Glx',target{~tmp}};
        end

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

        fGM  = MRS_struct.out.(vox{kk}).tissue.fGM(ii);
        fWM  = MRS_struct.out.(vox{kk}).tissue.fWM(ii);
        fCSF = MRS_struct.out.(vox{kk}).tissue.fCSF(ii);

        % Gasparovic et al. method (RAEE)
        % Calculate molal fractions from volume fractions (equivalent to eqs. 5-7 in Gasparovic et al., 2006)
        molal_fGM  = (fGM * concW_GM) / (fGM * concW_GM + fWM * concW_WM + fCSF * concW_CSF);
        molal_fWM  = (fWM * concW_WM) / (fGM * concW_GM + fWM * concW_WM + fCSF * concW_CSF);
        molal_fCSF = (fCSF * concW_CSF) / (fGM * concW_GM + fWM * concW_WM + fCSF * concW_CSF);

        for jj = 1:length(target)

            switch target{jj}

                case 'GABA'
                    EditingEfficiency = 0.5; % For TE = 68 ms
                    T1_Metab  = 1.31;  % Puts et al. 2013 (JMRI)
                    T2_Metab  = 0.088; % Edden et al. 2012 (JMRI)
                    N_H_Metab = 2;
                    MM  = 0.45; % MM correction: fraction of GABA in GABA+ peak. (In TrypDep, 30 subjects: 55% of GABA+ was MM)
                                % This fraction is platform- and implementation-dependent, based on length and
                                % shape of editing pulses and ifis Henry method
                    cWM = 1; % relative intrinsic concentration of GABA in pure WM
                    cGM = 2; % relative intrinsic concentration of GABA in pure GM

                case 'Glx'
                    EditingEfficiency = 0.4; % determined by FID-A simulations (for TE = 68 ms)
                    T1_Metab  = 1.23; % Posse et al. 2007 (MRM)
                    T2_Metab  = 0.18; % Ganji et al. 2012 (NMR Biomed)
                    N_H_Metab = 1;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of Glx in pure WM
                    cGM = 2; % relative intrinsic concentration of Glx in pure GM

                case 'GSH'
                    EditingEfficiency = 0.74; % At 3T based on Quantification of Glutathione in the Human Brain by MR Spectroscopy at 3 Tesla:
                                              % Comparison of PRESS and MEGA-PRESS
                                              % Faezeh Sanaei Nezhad etal. DOI 10.1002/mrm.26532, 2016
                    T1_Metab  = 0.40; % At 3T based on Doubly selective multiple quantum chemical shift imaging and
                                      % T1 relaxation time measurement of glutathione (GSH) in the human brain in vivo
                                      % In-Young Choi et al. NMR Biomed. 2013; 26: 28-34
                    T2_Metab  = 0.12; % At 3T based on the ISMRM abstract
                                      % T2 relaxation times of 18 brain metabolites determined in 83 healthy volunteers in vivo
                                      % Milan Scheidegger et al. Proc. Intl. Soc. Mag. Reson. Med. 22 (2014)
                    N_H_Metab = 2;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of GSH in pure WM
                    cGM = 1; % relative intrinsic concentration of GSH in pure GM

                case 'Lac'
                    EditingEfficiency = 0.94; % determined by FID-A simulations (for TE = 140 ms)
                    T1_Metab  = 1.50; % Wijnen et al. 2015 (NMR Biomed)
                    T2_Metab  = 0.24; % Madan et al. 2015 (MRM) (NB: this was estimated in brain tumors)
                    N_H_Metab = 3;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of Lac in pure WM
                    cGM = 1; % relative intrinsic concentration of Lac in pure GM

                case 'EtOH'
                    EditingEfficiency = 0.5; % assuming same as GABA for now
                    T1_Metab  = 1.31;  % assuming same as GABA
                    T2_Metab  = 0.088; % assuming same as GABA
                    N_H_Metab = 3;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of EtOH in pure WM
                    cGM = 1; % relative intrinsic concentration of EtOH in pure GM

                case 'Cr' % 3 ppm moiety
                    EditingEfficiency = 1; % not edited, so 1
                    T1_Metab  = (1.46 + 1.24)/2; % Mlynárik et al. 2001 (NMR in Biomed)
                    T2_Metab  = (166 + 144 + 148)/3/1e3; % Wyss et al. 2018 (MRM)
                    N_H_Metab = 3;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of Cr in pure WM
                    cGM = 1.5; % relative intrinsic concentration of Cr in pure GM

                case 'Cho' % 3.2 ppm moiety
                    EditingEfficiency = 1; % not edited, so 1
                    T1_Metab  = (1.30 + 1.08)/2; % Mlynárik et al. 2001 (NMR in Biomed)
                    T2_Metab  = (218 + 222 + 274)/3/1e3; % Wyss et al. 2018 (MRM)
                    N_H_Metab = 9;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of Cho in pure WM
                    cGM = 1; % relative intrinsic concentration of Cho in pure GM

                case 'NAA' % 2 ppm moiety
                    EditingEfficiency = 1; % not edited, so 1
                    T1_Metab  = (1.47 + 1.35)/2; % Mlynárik et al. 2001 (NMR in Biomed)
                    T2_Metab  = (343 + 263 + 253)/3/1e3; % Wyss et al. 2018 (MRM)
                    N_H_Metab = 3;
                    MM  = 1;
                    cWM = 1; % relative intrinsic concentration of NAA in pure WM
                    cGM = 1.5; % relative intrinsic concentration of NAA in pure GM

            end

            % Gasparovic et al. method (RAEE)
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_TissCorr(ii) = ...
                (MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) ./ MRS_struct.out.(vox{kk}).water.Area(ii)) .* ...
                (N_H_Water ./ N_H_Metab) .* MM ./ EditingEfficiency .* molal_concW .* ...
                (molal_fGM  .* (1 - exp(-TR_water./T1w_GM)) .* exp(-TE_water./T2w_GM) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab)) + ...
                molal_fWM  .* (1 - exp(-TR_water./T1w_WM)) .* exp(-TE_water./T2w_WM) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab)) + ...
                molal_fCSF .* (1 - exp(-TR_water./T1w_CSF)) .* exp(-TE_water./T2w_CSF) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab))) ./ ...
                (1 - molal_fCSF);

            % Alpha correction (Harris et al., 2015, JMRI)
            alpha = cWM ./ cGM;
            GrpAvgNorm = (meanfGM + alpha .* meanfWM) ./ ((fGM + alpha .* fWM) .* (meanfGM + meanfWM));
            ConcIU_TissCorr_Harris = ...
                (MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) ./ MRS_struct.out.(vox{kk}).water.Area(ii)) .* ...
                (N_H_Water ./ N_H_Metab) .* MM ./ EditingEfficiency .* ...
                (fGM .* concW_GM .* (1 - exp(-TR_water./T1w_GM)) .* exp(-TE_water./T2w_GM) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab)) + ...
                fWM .* concW_WM .* (1 - exp(-TR_water./T1w_WM)) .* exp(-TE_water./T2w_WM) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab)) + ...
                fCSF .* concW_CSF .* (1 - exp(-TR_water./T1w_CSF)) .* exp(-TE_water./T2w_CSF) ./ ((1 - exp(-TR./T1_Metab)) .* exp(-TE./T2_Metab)));
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr(ii)         = ConcIU_TissCorr_Harris ./ (fGM + alpha .* fWM);
            MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr_GrpNorm(ii) = ConcIU_TissCorr_Harris .* GrpAvgNorm;
            MRS_struct.out.(vox{kk}).(target{jj}).alpha = alpha;

        end

        target = target(1:end-3); % Remove Cr, Cho, and NAA from next steps
        
        % Build output figure
        if ishandle(105)
            clf(105);
        end
        if MRS_struct.p.hide
            h = figure('Visible', 'off');
        else
            h = figure(105);
        end
        % Open figure in center of screen
        scr_sz = get(0,'ScreenSize');
        fig_w = 1000;
        fig_h = 707;
        set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
        set(h,'Color',[1 1 1]);
        figTitle = 'GannetQuantify Output';
        set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');

        % Segmented voxel montage
        img_montage = MRS_struct.mask.(vox{kk}).img_montage{ii};
        img_montage = [img_montage(:,1:size(img_montage,2)/2,:); img_montage(:,size(img_montage,2)/2+1:end,:)];

        ha = subplot(2,2,1);
        imagesc(img_montage);
        colormap('gray');
        img = MRS_struct.mask.(vox{kk}).img{ii}(:);
        caxis([0 mean(img(img > 0.01)) + 3*std(img(img > 0.01))]); %#ok<*CAXIS>
        axis equal tight off;
        pos = get(ha,'pos');
        s = 0.04;
        set(ha,'pos',[pos(1)-s, pos(2)-s-0.02, pos(3)+2*s, pos(4)+2*s]);

        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,tmp,tmp2] = fileparts(MRS_struct.metabfile{1,ii});
        end
        fname = [tmp tmp2];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        [~,tmp3,tmp4] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp3 tmp4];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        title(sprintf(['Voxel from ' fname ' on ' T1image]), 'Interpreter', 'none');

        % Post-alignment spectra + model fits
        subplot(2,2,3);
        PlotPrePostAlign2(MRS_struct, vox, ii);

        % Output results
        subplot(2,2,2);
        axis off;

        tmp = strcmp(target,'GABAGlx');
        if any(tmp)
            if MRS_struct.p.HERMES
                target = {'GABA','Glx',target{~tmp}};
            else
                target = {'GABA','Glx'};
            end
        end

        text_pos = 1;

        tmp1 = 'Filename: ';
        if strcmp(MRS_struct.p.vendor,'Siemens_rda')
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{1,ii*2-1});
        else
            [~,tmp2,tmp3] = fileparts(MRS_struct.metabfile{1,ii});
        end
        fname = [tmp2 tmp3];
        if length(fname) > 30
            fname = [fname(1:12) '...' fname(end-11:end)];
        end
        text(0.4, text_pos-0.1, tmp1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        if MRS_struct.p.join
            text(0.425, text_pos-0.1, [fname ' (+ ' num2str(MRS_struct.p.numFilesPerScan - 1) ' more)'], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
        else
            text(0.425, text_pos-0.1, fname, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
        end

        tmp1 = 'Anatomical image: ';
        [~,tmp2,tmp3] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
        T1image = [tmp2 tmp3];
        if length(T1image) > 30
            T1image = [T1image(1:12) '...' T1image(end-11:end)];
        end
        text(0.4, text_pos-0.2, tmp1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        text(0.425, text_pos-0.2, T1image, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');

        for jj = 1:length(target)

            switch target{jj}
                case 'GABA'
                    tmp2 = 'GABA+';
                case 'Lac'
                    tmp2 = 'Lac+MM';
                case {'Glx','GSH','EtOH'}
                    tmp2 = target{jj};
            end

            shift = 0;
            for ll = 1:3

                text_pos = 0.7;
                if ll == 1
                    tmp1 = 'Relaxation-, tissue-corrected (Gasparovic et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_TissCorr(ii));
                    text(0, text_pos, tmp1, 'Units', 'normalized', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                elseif ll == 2
                    text_pos = text_pos - 0.2 - shift;
                    tmp1 = 'Relaxation-, tissue-, alpha-corrected (Harris et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr(ii));
                    text(0, text_pos, tmp1, 'Units', 'normalized', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                elseif ll == 3
                    text_pos = text_pos - 0.4 - shift;
                    tmp1a = 'Relaxation-, tissue-, alpha-corrected; group-average-normalized';
                    tmp1b = '(Harris et al. method)';
                    tmp3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU_AlphaTissCorr_GrpNorm(ii));
                    text(0, text_pos, tmp1a, 'Units', 'normalized', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                    text(0, text_pos - 0.1, tmp1b, 'Units', 'normalized', 'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
                end
                if ll == 3
                    text_pos = text_pos - 0.1*(jj+1);
                else
                    text_pos = text_pos - 0.1*jj;
                end
                if ll == 1
                    text(0.4, text_pos, [tmp2 '/Water: '], 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                else
                    text(0.4, text_pos, [tmp2 '/Water (\alpha = ' num2str(MRS_struct.out.(vox{kk}).(target{jj}).alpha) '): '], ...
                        'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
                end
                text(0.425, text_pos, tmp3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10);

                if MRS_struct.p.HERMES
                    shift = shift + 0.1*(numel(target)-1);
                else
                    shift = shift + 0.1;
                end

            end

        end

        text(0.4, text_pos - 0.15, 'QuantifyVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10, 'HorizontalAlignment', 'right');
        text(0.425, text_pos - 0.15, MRS_struct.version.quantify, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10);

        % Save output as PDF
        run_count = SavePDF(h, MRS_struct, ii, 1, kk, vox, mfilename, run_count);

    end

    if MRS_struct.p.mat % save MRS_struct as mat file
        mat_name = fullfile(pwd, ['MRS_struct_' vox{kk} '.mat']);
        if exist(mat_name, 'file')
            fprintf('\nUpdating results in %s\n', ['MRS_struct_' vox{kk} '.mat...']);
        else
            fprintf('\nSaving results to %s\n', ['MRS_struct_' vox{kk} '.mat...']);
        end
        save(mat_name, 'MRS_struct', '-v7.3');
    end

    if MRS_struct.p.csv % export MRS_struct fields into csv file
        MRS_struct = ExportToCSV(MRS_struct, vox{kk}, 'quantify');
    end

end

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end



