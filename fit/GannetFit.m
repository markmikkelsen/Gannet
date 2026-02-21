function MRS_struct = GannetFit(MRS_struct, varargin)
% Signal fitting in the frequency domain using nonlinear least-squares optimization

if nargin == 0
    fprintf('\n');
    error('MATLAB:minrhs', 'Not enough input arguments.');
end

if ~isstruct(MRS_struct)
    fprintf('\n');
    error('The first input argument must be a structure, but received %s.', class(MRS_struct));
end

MRS_struct.info.datetime.fit = datetime('now');
MRS_struct.info.version.fit = '260217';

if ~isMATLABReleaseOlderThan("R2025a") && MRS_struct.p.append
    font_size_adj  = 2.75;
    font_size_adj2 = 1.5;
else
    font_size_adj  = 0;
    font_size_adj2 = 0;
end

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end

if MRS_struct.p.phantom
    error('The loaded data are phantom data. Use GannetFitPhantom instead of GannetFit.');
end

if nargin > 1
    % varargin = Optional arguments if user wants to specify a target
    % metabolite, overwriting the parameter set in GannetPreInitialise.m
    switch varargin{1}
        case 'GABA'
            MRS_struct.p.target = {'GABA'};
        case 'Glx'
            MRS_struct.p.target = {'Glx'};
        case 'GABAGlx'
            MRS_struct.p.target = {'GABAGlx'};
        case 'GSH'
            MRS_struct.p.target = {'GSH'};
        case 'Lac'
            MRS_struct.p.target = {'Lac'};
        case 'EtOH'
            MRS_struct.p.target = {'EtOH'};
    end
end

target = MRS_struct.p.target;
freq = MRS_struct.spec.freq;

lsqopts = optimset('lsqcurvefit');
lsqopts = optimset(lsqopts, 'MaxIter', 800, 'TolX', 1e-4, 'TolFun', 1e-4, 'Display', 'off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts, 'MaxIter', 800, 'TolX', 1e-6, 'TolFun', 1e-6, 'FunValCheck', 'off');
lsqnlinopts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxFunctionEvaluations', 1e4, ...
                           'MaxIterations', 1e3, 'FunctionTolerance', 1e-7);

warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:IterationLimitExceeded');
warning('off','MATLAB:rankDeficientMatrix');

% Loop over voxels if PRIAM
for kk = 1:length(vox)

    if strcmp(MRS_struct.p.reference, 'H2O')
        WaterData = MRS_struct.spec.(vox{kk}).water;
    end

    run_count = 0;

    % Loop over edited spectra if HERMES
    for jj = 1:length(target)

        if strcmp(target{jj}, 'GABAGlx')
            t = 'GABA+Glx';
        else
            t = target{jj};
        end

        if jj == 1
            fprintf('\nFitting %s in...\n', t);
        else
            fprintf('Fitting %s in...\n', t);
        end

        DIFF = MRS_struct.spec.(vox{kk}).(target{jj}).diff;
        SUM  = MRS_struct.spec.(vox{kk}).(target{jj}).sum;
        if MRS_struct.p.HERMES
            OFF = MRS_struct.spec.(vox{kk}).(target{jj}).off_off;
        else
            OFF = MRS_struct.spec.(vox{kk}).(target{jj}).off;
        end

        error_report = cell(1);
        catch_ind    = 1;

        for ii = 1:MRS_struct.p.numScans

            [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii});
            fprintf('%s...\n', [name ext]);

            try % pass to next dataset if errors occur

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1.  Metabolite Fitting
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                switch target{jj}
                    case 'GABAGlx'
                        fitFun = @FitGABAGlx;
                    case 'GABA'
                        fitFun = @FitGABA;
                    case 'Glx'
                        fitFun = @FitGlx;
                    case 'GSH'
                        fitFun = @FitGSH;
                    case 'Lac'
                        fitFun = @FitLac;
                    case 'EtOH'
                        fitFun = @FitEtOH;
                    otherwise
                        error('Metabolite ''%s'' not recognized.', target{jj});
                end

                % Baseline modeling
                window_size = floor(1./MRS_struct.p.SpecResNominal(ii)); % 1-Hz window size
                lambda_DIFF = 10.^(floor(log(length(DIFF(ii,:)))) + 2); % log(lambda) for Whittaker smoother is roughly proportional to log(N_datapoints)
                lambda_SUM  = 10.^(floor(log(length(DIFF(ii,:)))) - 1);

                DIFF_tmp = real(DIFF(ii,:));
                base_mask = BaselineRecognition(DIFF_tmp, freq, window_size);
                baseline.DIFF = BaselineSmoothing(ii*2-1, freq, DIFF_tmp, base_mask, lambda_DIFF).';

                SUM_tmp = real(SUM(ii,:));
                base_mask = BaselineRecognition(SUM_tmp, freq, window_size);
                baseline.SUM = BaselineSmoothing(ii*2, freq, SUM_tmp, base_mask, lambda_SUM).';

                h_tmp = figure('Visible','off');
                % h_tmp = figure(333);
                clf(h_tmp);
                set(h_tmp,'Units','Normalized','OuterPosition',[0 0 0.5 1]);
                tiledlayout(2,1);

                if strcmp(target{jj}, 'GABAGlx')
                    xlims = [2.075 4.5];
                else
                    xlims = [0.5 4.5];
                end

                nexttile;
                hold on;
                plot(freq, DIFF_tmp, 'k', 'LineWidth', 1);
                plot(freq, baseline.DIFF, 'r', 'LineWidth', 1);
                hold off;
                set(gca,'XDir','reverse','TickDir','out','XLim',xlims,'FontSize',14);
                ylims = ylim;
                set(gca,'XLim',[-4 10],'YLim',ylims);
                title([target{jj} '-edited'],'FontSize',18);
                
                nexttile;
                hold on;
                plot(freq, SUM_tmp, 'k', 'LineWidth', 1);
                plot(freq, baseline.SUM, 'r', 'LineWidth', 1);
                hold off;
                xlabel('ppm','FontSize',16,'FontWeight','bold');
                set(gca,'XDir','reverse','TickDir','out','XLim',[0.5 4.5],'FontSize',14);
                ylims = ylim;
                set(gca,'XLim',[-2 8],'YLim',ylims);
                title('SUM','FontSize',18);
                
                legend({'data','baseline'},'Box','off','Location','best','FontSize',14);
                drawnow;

                out_dir = fullfile(pwd, 'Gannet_model_output');
                if ~exist(out_dir, 'dir')
                    mkdir(out_dir);
                end

                [~,fname,ext] = fileparts(MRS_struct.metabfile{ii});
                if strcmpi(ext, '.gz')
                    fname(end-3:end) = [];
                end
                exportgraphics(h_tmp, fullfile(out_dir, [fname '_' target{jj} '_baseline_model.png']), "Resolution", 300);
                close(h_tmp);

                % Fit metabolite signal model
                [MRS_struct, modelFit] = fitFun(MRS_struct, freq, DIFF, vox, target, ii, jj, kk, baseline.DIFF, lsqopts, nlinopts, lsqnlinopts);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1a. Initialize the Output Figure
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ishandle(102)
                    clf(102);
                end
                if MRS_struct.p.hide
                    h = figure('Visible', 'off');
                else
                    h = figure(102);
                end
                if ~isMATLABReleaseOlderThan("R2025a")
                    h.Theme = 'light';
                end
                % Open figure in center of screen
                scr_sz = get(0,'ScreenSize');
                fig_w  = 1000;
                fig_h  = 707;
                set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
                set(h,'Color',[1 1 1]);
                figTitle = 'GannetFit Output';
                set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
                
                % Spectra plot
                subplot(2,2,1);
                metabMin  = min(real(DIFF(ii,modelFit.plotBounds)));
                metabMax  = max(real(DIFF(ii,modelFit.plotBounds)));
                resMax    = max(modelFit.residPlot);
                residPlot = modelFit.residPlot + metabMin - resMax;
                
                switch target{jj}
                    case 'GABA'
                        residPlot = residPlot + metabMin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(modelFit.weightRange) = NaN;
                        hold on;
                        plot(freq(modelFit.plotBounds), real(DIFF(ii,modelFit.plotBounds)), 'b', ...
                            freq(modelFit.freqBounds), GaussModel(modelFit.full.modelParam, freq(modelFit.freqBounds)), 'r', ...
                            freq(modelFit.freqBounds), residPlot2, 'k');
                        plot(freq(modelFit.freqBounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                        hold off;
                        set(gca, 'XLim', [2.6 3.6], 'FontSize', 10 - font_size_adj);
                        
                    case 'GSH'
                        residPlot = residPlot + metabMin - max(residPlot);
                        residPlot2 = residPlot;
                        GSHgaussModel = @EightGaussModel_noBaseline;
                        hold on;
                        plot(freq(modelFit.plotBounds), real(DIFF(ii,modelFit.plotBounds)), 'b' , ...
                            freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.full, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds), 'r', ...
                            freq(modelFit.freqBounds), residPlot2, 'k');
                        if MRS_struct.p.show_fits
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss1, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss2, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss3, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss4, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss5, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss6, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss7, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GSHgaussModel(modelFit.modelParam.Gauss8, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                        end
                        hold off;
                        set(gca, 'XLim', [1.8 4.2], 'FontSize', 10 - font_size_adj);

                    case 'Lac'
                        hold on;
                        p1 = plot(freq(plotBounds), real(DIFF(ii,plotBounds)), 'k');
                        p2 = plot(freq(freqBounds), LacModel(BHBmodelParam,freq(freqBounds)), 'Color', [31 120 180]/255, 'LineWidth', 1);
                        p3 = plot(freq(freqBounds), LacModel(LacModelParam,freq(freqBounds)), 'Color', [51 160 44]/255, 'LineWidth', 1);
                        p4 = plot(freq(freqBounds), LacModel(LacPlusModelParam,freq(freqBounds)), 'r', 'LineWidth', 1);
                        p5 = plot(freq(freqBounds), residPlot, 'k');
                        hold off;
                        set(gca, 'XLim', [0.7 1.9], 'XTick', 0:0.25:10, 'FontSize', 10 - font_size_adj);
                        
                    case 'Glx'
                        plot(freq(plotBounds), real(DIFF(ii,plotBounds)), 'b', ...
                            freq(freqBounds), DoubleGaussModel(GaussModelParam,freq(freqBounds)), 'r', ...
                            freq(freqBounds), residPlot, 'k');
                        set(gca, 'XLim', [3.4 4.2], 'FontSize', 10 - font_size_adj);
                        
                    case 'GABAGlx'
                        residPlot  = residPlot + metabMin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(modelFit.weightRange) = NaN;
                        hold on;
                        plot(freq(modelFit.plotBounds), real(DIFF(ii,modelFit.plotBounds)), 'b', ...
                            freq(modelFit.plotBounds), GABAGlxModel_noBaseline(modelFit.modelParam.full, freq(modelFit.plotBounds)) + ...
                                baseline.DIFF(modelFit.plotBounds), 'r', ...
                            freq(modelFit.freqBounds), residPlot2, 'k');
                        % Plot weighted portion of residuals in different color
                        if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor,{'Philips','Philips_data','Philips_raw'}))
                            plot(freq(modelFit.freqBounds(modelFit.ChoRange)), residPlot(modelFit.ChoRange), 'Color', [255 160 64]/255);
                            plot(freq(modelFit.freqBounds(modelFit.GlxDownfieldRange)), residPlot(modelFit.GlxDownfieldRange), 'Color', [255 160 64]/255);
                        else
                            plot(freq(modelFit.freqBounds(modelFit.weightRange)), residPlot(modelFit.weightRange), 'Color', [255 160 64]/255);
                        end
                        if MRS_struct.p.show_fits
                            plot(freq(modelFit.freqBounds), GABAGlxModel_noBaseline(modelFit.modelParam.Gauss1, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GABAGlxModel_noBaseline(modelFit.modelParam.Gauss2, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                            plot(freq(modelFit.freqBounds), GABAGlxModel_noBaseline(modelFit.modelParam.Gauss3, freq(modelFit.freqBounds)) + ...
                                baseline.DIFF(modelFit.freqBounds));
                        end
                        hold off;
                        set(gca, 'XLim', [2.6 4.2], 'XTick', 0:0.2:10, 'FontSize', 10 - font_size_adj);

                    case 'EtOH'
                        residPlot = residPlot + metabMin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(weightRange) = NaN;
                        hold on;
                        plot(freq(plotBounds), real(DIFF(ii,plotBounds)), 'b', ...
                            freq(freqBounds), EtOHModel(LorentzModelParam,freq(freqBounds)), 'r', ...
                            freq(freqBounds), residPlot2, 'k');
                        plot(freq(freqBounds(LacRange)), residPlot(LacRange), 'Color', [255 160 64]/255);
                        hold off;
                        set(gca, 'XLim', [0.4 1.9], 'FontSize', 10 - font_size_adj);
                end
                
                % From here on is cosmetic - adding labels etc.
                switch target{jj}
                    case 'GABA'
                        text(3, metabMax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelBounds = freq <= 2.4 & freq >= 2;
                        tailTop = max(real(DIFF(ii,labelBounds)));
                        tailBottom = min(real(DIFF(ii,labelBounds)));
                        text(2.8, min(residPlot), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(2.8, tailTop + metabMax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.8, tailBottom - metabMax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        if length(MRS_struct.p.target) ~= 3
                            text(3.2, min(residPlot) - 0.5*abs(max(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                                'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        end
                        
                    case 'GSH'
                        text(2.95, metabMax+metabMax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelBounds = freq <= 2.4 & freq >= 1.75;
                        tailTop = max(real(DIFF(ii,labelBounds)));
                        tailBottom = min(real(DIFF(ii,labelBounds)));
                        text(2.25, min(residPlot), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(2.25, tailTop + metabMax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.45, tailBottom - 20*metabMax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            text(3.3, max(residPlot) + 0.5*abs(max(residPlot) - min(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                                'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
                        end
                        
                    case 'Lac'
                        % text(1.27, metabmax/3.5, 'Lac+', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        % labelbounds = freq <= 0.8 & freq >= 0.42;
                        % tailtop = max(real(DIFF(ii,labelbounds)));
                        % tailbottom = min(real(DIFF(ii,labelbounds)));
                        % text(0.44, mean(real(DIFF(ii,labelbounds))) + 4*std(real(DIFF(ii,labelbounds))), 'data', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
                        % text(0.44, mean(residPlot) - 4*std(residPlot), 'residual', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
                        % text(0.85, tailbottom, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        legend([p1, p2, p3, p4, p5], {'data','BHB+','Lac+','model','residual'}, 'box', 'off', 'Location', 'northeast', 'FontSize', 9 - font_size_adj);
                        
                    case 'Glx'
                        text(3.8, metabMax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelBounds = freq <= 3.6 & freq >= 3.4;
                        tailTop = max(real(DIFF(ii,labelBounds)));
                        tailBottom = min(real(DIFF(ii,labelBounds)));
                        text(3.5, min(residPlot),'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(3.5, tailTop + metabMax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(3.5, tailBottom - metabMax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        
                    case 'GABAGlx'
                        text(3, metabMax/4, 'GABA', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        text(3.755, metabMax/4, 'Glx', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        text(2.8, min(residPlot), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        labelBounds = freq <= 2.8 & freq >= 2.7;
                        tailTop = max(real(DIFF(ii,labelBounds)));
                        tailBottom = min(real(DIFF(ii,labelBounds)));
                        text(2.8, tailTop + metabMax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.8, tailBottom - metabMax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        text(3.2, min(residPlot) - 0.5*abs(max(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                            'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        
                    case 'EtOH'
                        text(1.45, metabMax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelBounds = freq <= 0.9 & freq >= 0.5;
                        tailTop = max(real(DIFF(ii,labelBounds)));
                        tailBottom = min(real(DIFF(ii,labelBounds)));
                        text(0.6, min(residPlot),'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(0.8, tailTop + metabMax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(0.8, tailBottom - metabMax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        text(1.4, max(residPlot2) + 0.5*abs(max(residPlot2)), 'down-weighted', 'Color', [255 160 64]/255, ...
                            'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                end
                
                title('Difference spectrum and model fit', 'FontSize', 11 - font_size_adj);
                xlabel('ppm', 'FontSize', 11 - font_size_adj);
                set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'FontSize', 10 - font_size_adj);
                set(get(gca,'YAxis'), 'Visible', 'off');
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   2.  Water Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if strcmp(MRS_struct.p.reference,'H2O')
                    [MRS_struct, hb] = FitWater(MRS_struct, freq, WaterData, OFF, vox, target, ii, jj, kk, lsqopts, nlinopts);
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   3.  Cr, NAA, and Glu Fits
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [MRS_struct, Cr_OFF, freqBoundsChoCr, freqBoundsCr, residCr] = ...
                    FitCrNAAGlu(MRS_struct, freq, OFF, SUM, vox, ii, kk, baseline.SUM, lsqopts, nlinopts, lsqnlinopts);


                % Root sum square fit errors and concentrations as metabolite ratios
                if strcmpi(target{jj},'GABAGlx')
                    MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii)  = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii)   = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).GABA.FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).Glx.FitError_NAA(ii)  = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).GABA.ConcCr(ii)       = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                    MRS_struct.out.(vox{kk}).GABA.ConcCho(ii)      = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                    MRS_struct.out.(vox{kk}).GABA.ConcNAA(ii)      = MRS_struct.out.(vox{kk}).GABA.Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
                    MRS_struct.out.(vox{kk}).Glx.ConcCr(ii)        = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                    MRS_struct.out.(vox{kk}).Glx.ConcCho(ii)       = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                    MRS_struct.out.(vox{kk}).Glx.ConcNAA(ii)       = MRS_struct.out.(vox{kk}).Glx.Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
                else
                    MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii)  = sqrt(MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).(target{jj}).FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii)       = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                    MRS_struct.out.(vox{kk}).(target{jj}).ConcCho(ii)      = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                    MRS_struct.out.(vox{kk}).(target{jj}).ConcNAA(ii)      = MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);
                end
                
                MRS_struct.out.(vox{kk}).Cho.FitError_Cr(ii)  = sqrt(MRS_struct.out.(vox{kk}).Cho.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).NAA.FitError_Cr(ii)  = sqrt(MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).Glu.FitError_Cr(ii)  = sqrt(MRS_struct.out.(vox{kk}).Glu.FitError(ii).^2 + MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2);
                MRS_struct.out.(vox{kk}).Glu.FitError_NAA(ii) = sqrt(MRS_struct.out.(vox{kk}).Glu.FitError(ii).^2 + MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2);

                MRS_struct.out.(vox{kk}).Cho.ConcCr(ii)  = MRS_struct.out.(vox{kk}).Cho.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).NAA.ConcCr(ii)  = MRS_struct.out.(vox{kk}).NAA.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).Glu.ConcCr(ii)  = MRS_struct.out.(vox{kk}).Glu.Area(ii) / MRS_struct.out.(vox{kk}).Cr.Area(ii);
                MRS_struct.out.(vox{kk}).Glu.ConcCho(ii) = MRS_struct.out.(vox{kk}).Glu.Area(ii) / MRS_struct.out.(vox{kk}).Cho.Area(ii);
                MRS_struct.out.(vox{kk}).Glu.ConcNAA(ii) = MRS_struct.out.(vox{kk}).Glu.Area(ii) / MRS_struct.out.(vox{kk}).NAA.Area(ii);

                if strcmp(MRS_struct.p.reference,'H2O')
                    MRS_struct.out.(vox{kk}).Cr.FitError_W(ii)  = sqrt(MRS_struct.out.(vox{kk}).Cr.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).Cho.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Cho.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).NAA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).NAA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                    MRS_struct.out.(vox{kk}).Glu.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glu.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);

                    MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Cr', ii);
                    MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Cho', ii);
                    MRS_struct = CalcIU(MRS_struct, vox{kk}, 'NAA', ii);
                    MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glu', ii);
                end

                % Reorder structure fields
                if strcmp(MRS_struct.p.reference,'H2O')
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'FitError_NAA', 'FitError_W', 'ConcCr', 'ConcCho', 'ConcNAA', 'ConcIU'};
                    if strcmpi(target{jj},'GABAGlx')
                        MRS_struct.out.(vox{kk}).GABA = orderfields(MRS_struct.out.(vox{kk}).GABA, fields);
                        MRS_struct.out.(vox{kk}).Glx  = orderfields(MRS_struct.out.(vox{kk}).Glx, fields);
                    else
                        MRS_struct.out.(vox{kk}).(target{jj}) = orderfields(MRS_struct.out.(vox{kk}).(target{jj}), fields);
                    end
                else
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'FitError_NAA', 'ConcCr', 'ConcCho', 'ConcNAA'};
                    if strcmpi(target{jj},'GABAGlx')
                        MRS_struct.out.(vox{kk}).GABA = orderfields(MRS_struct.out.(vox{kk}).GABA, fields);
                        MRS_struct.out.(vox{kk}).Glx  = orderfields(MRS_struct.out.(vox{kk}).Glx, fields);
                    else
                        MRS_struct.out.(vox{kk}).(target{jj}) = orderfields(MRS_struct.out.(vox{kk}).(target{jj}), fields);
                    end
                end

                if strcmp(MRS_struct.p.reference,'H2O')
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_W', 'ConcIU'};
                else
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError'};
                end
                MRS_struct.out.(vox{kk}).Cr  = orderfields(MRS_struct.out.(vox{kk}).Cr, fields);

                if strcmp(MRS_struct.p.reference,'H2O')
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'FitError_W', 'ConcCr', 'ConcIU'};
                else
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'ConcCr'};
                end
                MRS_struct.out.(vox{kk}).Cho = orderfields(MRS_struct.out.(vox{kk}).Cho, fields);
                MRS_struct.out.(vox{kk}).NAA = orderfields(MRS_struct.out.(vox{kk}).NAA, fields);

                if strcmp(MRS_struct.p.reference,'H2O')
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'FitError_NAA', 'FitError_W', 'ConcCr', 'ConcCho', 'ConcNAA', 'ConcIU'};
                else
                    fields = {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError', 'FitError_Cr', 'FitError_NAA', 'ConcCr', 'ConcCho', 'ConcNAA'};
                end
                MRS_struct.out.(vox{kk}).Glu = orderfields(MRS_struct.out.(vox{kk}).Glu, fields);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   6. Build GannetFit Output
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                CrMin = min(real(Cr_OFF(freqBoundsCr)));
                CrMax = max(real(Cr_OFF(freqBoundsCr)));
                resMaxCr = max(residCr);
                residCr = residCr + CrMin - resMaxCr;
                if strcmp(MRS_struct.p.reference,'H2O')
                    hb_pos = get(hb, 'Position');
                    hd = axes(hb.Parent, 'Position', [1.3*hb_pos(1)+hb_pos(3)-0.175, 1.001*hb_pos(2)+hb_pos(4)-0.175 0.1125 0.1575]);
                    plot(freq, real(Cr_OFF), 'b', ...
                        freq(freqBoundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqBoundsChoCr))), 'r', ...
                        freq(freqBoundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqBoundsChoCr))), 'r', ...
                        freq(freqBoundsCr), residCr, 'k');
                    text(2.94, CrMax*0.75, 'Creatine', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj2);
                    set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XLim', [2.6 3.6], 'XTick', 2.6:0.2:3.6, 'FontSize', 6 - font_size_adj2);
                    xlabel('ppm', 'FontSize', 11 - font_size_adj);
                    set(get(gca,'YAxis'), 'Visible', 'off');
                    hd_pos = get(hd, 'Position');
                    set(hd, 'Position', [hb_pos(1)+hb_pos(3)-hd_pos(3) hd_pos(2:4)]);
                else
                    hb = subplot(2,2,3);
                    plot(freq, real(Cr_OFF), 'b', ...
                        freq(freqBoundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqBoundsChoCr))), 'r', ...
                        freq(freqBoundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqBoundsChoCr))), 'r', ...
                        freq(freqBoundsCr), residCr, 'k');
                    set(gca,'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XTick', 2.6:0.2:3.6, 'FontSize', 10 - font_size_adj);
                    set(get(gca,'YAxis'),'Visible','off');
                    xlim([2.6 3.6]);
                    Crlabelbounds = freq(freqBoundsCr) <= 3.12 & freq(freqBoundsCr) >= 2.72;
                    text(3, max(residCr(Crlabelbounds)) + 0.05*CrMax, 'residual', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                    text(2.7, 0.1*CrMax, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                    text(2.7, 0.01*CrMax, 'model','Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                    text(2.94, CrMax*0.75, 'Creatine', 'FontSize', 10 - font_size_adj);
                    xlabel('ppm', 'FontSize', 11 - font_size_adj);
                    title('Reference signal', 'FontSize', 11 - font_size_adj);
                end
                
                if any(strcmp('mask',fieldnames(MRS_struct)))
                    hc = subplot(2,2,2);
                    set(hc,'Position',[0.52 0.52 0.42 0.42]) % move the axes slightly
                    size_max = size(MRS_struct.mask.img{ii},1);
                    imagesc(MRS_struct.mask.img{ii}(:,size_max+(1:size_max)));
                    colormap('gray');
                    caxis([0 1]); %#ok<*CAXIS> 
                    axis equal tight off;
                    subplot(2,2,4,'replace');
                else
                    subplot(3,2,[2 4]);
                    axis off;
                end
                
                % Cleaner text alignment; move GABA/Glx to separate lines
                text_pos = 0.95; % A variable to determine y-position of text on printout on figure
                shift = 0.06;
                
                % 1. Filename
                if strcmp(MRS_struct.p.vendor,'Siemens_rda')
                    [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii*2-1});
                else
                    [~,name,ext] = fileparts(MRS_struct.metabfile{1,ii});
                end
                fname = [name ext];
                if length(fname) > 30
                    fname = sprintf([fname(1:floor((end-1)/2)) '...\n     ' fname(ceil(end/2):end)]);
                    shift2 = 0.02;
                else
                    shift2 = 0;
                end
                text(0.4, text_pos, 'Filename: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                if MRS_struct.p.join
                    text(0.425, text_pos+shift2, [fname ' (+ ' num2str(MRS_struct.p.numFilesPerScan - 1) ' more)'], ...
                        'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'Interpreter', 'none');
                else
                    text(0.425, text_pos+shift2, fname, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'Interpreter', 'none');
                end
                
                % 2a. Area
                text(0.4, text_pos-shift, 'Area ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                
                switch target{jj}
                    case 'GABA'
                        str1 = 'GABA: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).GABA.Area(ii));
                    case 'Glx'
                        str1 = 'Glx: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Glx.Area(ii));
                    case 'GABAGlx'
                        str1 = 'GABA: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).GABA.Area(ii));
                        str3 = 'Glx: ';
                        str4 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Glx.Area(ii));
                    case 'GSH'
                        str1 = 'GSH: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
                    case 'Lac'
                        str1 = 'Lac+: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
                    case 'EtOH'
                        str1 = 'EtOH: ';
                        str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).(target{jj}).Area(ii));
                end
                
                text(0.4, text_pos-2*shift, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                text(0.425, text_pos-2*shift, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                
                if strcmp(target{jj},'GABAGlx')
                    text(0.4, text_pos-3*shift, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-3*shift, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    n = 0;
                else
                    n = shift;
                end
                
                if strcmp(MRS_struct.p.reference,'H2O')
                    
                    % 2b. Area (Water / Cr)
                    str1 = sprintf('%.3g', MRS_struct.out.(vox{kk}).water.Area(ii));
                    str2 = sprintf('%.3g', MRS_struct.out.(vox{kk}).Cr.Area(ii));
                    
                    text(0.4, text_pos-4*shift+n, 'Water: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-4*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    text(0.4, text_pos-5*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-5*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                    % 3. Linewidth
                    text(0.4, text_pos-6*shift+n, 'Linewidth ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    str1 = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).water.FWHM(ii));
                    str2 = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).Cr.FWHM(ii));
                    text(0.4, text_pos-7*shift+n, 'Water: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-7*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    text(0.4, text_pos-8*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-8*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                                        
                    % 4. SNR
                    text(0.4, text_pos-9*shift+n, 'SNR ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    str1 = sprintf('%.0f', MRS_struct.out.(vox{kk}).water.SNR(ii));
                    str2 = sprintf('%.0f', MRS_struct.out.(vox{kk}).Cr.SNR(ii));
                    text(0.4, text_pos-10*shift+n, 'Water: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-10*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    text(0.4, text_pos-11*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-11*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                    % 5a. Fit Error
                    text(0.4, text_pos-12*shift+n, 'Fit error ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    switch target{jj}
                        
                        case {'GABA','Glx','GSH','Lac','EtOH'}
                            
                            % 5b. Fit Error
                            if strcmpi(target{jj},'GABA')
                                str1 = 'GABA+,Water: ';
                                str2 = 'GABA+,Cr: ';
                            elseif strcmpi(target{jj},'Lac')
                                str1 = 'Lac+,Water: ';
                                str2 = 'Lac+,Cr: ';
                            else
                                str1 = sprintf('%s,Water: ', target{jj});
                                str2 = sprintf('%s,Cr: ', target{jj});
                            end
                            str3 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_W(ii));
                            str4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii));
                            
                            text(0.4, text_pos-13*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-13*shift+n, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-14*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-14*shift+n, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            % 6. Quantification
                            text(0.4, text_pos-15*shift+n, 'Quantification ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                            
                            if strcmpi(target{jj},'GABA')
                                str1 = 'GABA+/Water: ';
                                str2 = 'GABA+/Cr: ';
                            elseif strcmpi(target{jj},'Lac')
                                str1 = 'Lac+/Water: ';
                                str2 = 'Lac+/Cr: ';
                            else
                                str1 = sprintf('%s/Water: ', target{jj});
                                str2 = sprintf('%s/Cr: ', target{jj});
                            end
                            str3 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).(target{jj}).ConcIU(ii));
                            str4 = sprintf('%.2f', MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii));
                            
                            text(0.4, text_pos-16*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-16*shift+n, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-17*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-17*shift+n, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            n = 5*shift;
                            
                        case 'GABAGlx'
                            
                            % 5b. Fit Error
                            str1 = 'GABA+,Water: ';
                            str2 = 'GABA+,Cr: ';
                            str3 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_W(ii));
                            str4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii));
                            str5 = 'Glx,Water: ';
                            str6 = 'Glx,Cr: ';
                            str7 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_W(ii));
                            str8 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii));
                            
                            text(0.4, text_pos-13*shift, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-13*shift, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-14*shift, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-14*shift, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-15*shift, str5, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-15*shift, str7, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-16*shift, str6, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-16*shift, str8, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            % 6. Quantification
                            text(0.4, text_pos-17*shift, 'Quantification ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                            
                            str1 = 'GABA+/Water: ';
                            str2 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).GABA.ConcIU(ii));
                            str3 = 'GABA+/Cr: ';
                            str4 = sprintf('%.2f', MRS_struct.out.(vox{kk}).GABA.ConcCr(ii));
                            str5 = 'Glx/Water: ';
                            str6 = sprintf('%.2f i.u.', MRS_struct.out.(vox{kk}).Glx.ConcIU(ii));
                            str7 = 'Glx/Cr: ';
                            str8 = sprintf('%.2f', MRS_struct.out.(vox{kk}).Glx.ConcCr(ii));
                            
                            text(0.4, text_pos-18*shift, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-18*shift, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-19*shift, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-19*shift, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-20*shift, str5, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-20*shift, str6, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-21*shift, str7, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-21*shift, str8, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            n = 0;
                            
                    end
                    
                    % 7. FitVer
                    text(0.4, text_pos-22.5*shift+n, 'FitVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-22.5*shift+n, MRS_struct.info.version.fit, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                else
                    
                    % 2. Area (Cr)
                    str = sprintf('%.3g', MRS_struct.out.(vox{kk}).Cr.Area(ii));
                    
                    text(0.4, text_pos-4*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-4*shift+n, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                    % 3. Linewidth
                    text(0.4, text_pos-5*shift+n, 'Linewidth ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    str = sprintf('%.2f Hz', MRS_struct.out.(vox{kk}).Cr.FWHM(ii));
                    text(0.4, text_pos-6*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-6*shift+n, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                    % 4. SNR
                    text(0.4, text_pos-7*shift+n, 'SNR ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    str = sprintf('%.0f', MRS_struct.out.(vox{kk}).Cr.SNR(ii));
                    text(0.4, text_pos-8*shift+n, 'Cr: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-8*shift+n, str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                                        
                    % 5a. Fit Error
                    text(0.4, text_pos-9*shift+n, 'Fit error ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                    
                    switch target{jj}
                        
                        case {'GABA','Glx','GSH','Lac','EtOH'}
                            
                            % 4b. Fit Error
                            if strcmpi(target{jj},'GABA')
                                str1 = 'GABA+,Cr: ';
                            elseif strcmpi(target{jj},'Lac')
                                str1 = 'Lac+,Cr: ';
                            else
                                str1 = sprintf('%s,Cr: ', target{jj});
                            end
                            str2 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).(target{jj}).FitError_Cr(ii));
                            
                            text(0.4, text_pos-10*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-10*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            % 5. Quantification
                            text(0.4, text_pos-11*shift+n, 'Quantification ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                            
                            if strcmpi(target{jj},'GABA')
                                str1 = 'GABA+/Cr: ';
                            elseif strcmpi(target{jj},'Lac')
                                str1 = 'Lac+/Cr: ';
                            else
                                str1 = sprintf('%s/Cr: ', target{jj});
                            end
                            str2 = sprintf('%.2f', MRS_struct.out.(vox{kk}).(target{jj}).ConcCr(ii));
                            
                            text(0.4, text_pos-12*shift+n, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-12*shift+n, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            n = 3*shift;
                            
                        case 'GABAGlx'
                            
                            % 4b. Fit Error
                            str1 = 'GABA+,Cr: ';
                            str2 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).GABA.FitError_Cr(ii));
                            str3 = 'Glx,Cr: ';
                            str4 = sprintf('%.2f%%', MRS_struct.out.(vox{kk}).Glx.FitError_Cr(ii));
                            
                            text(0.4, text_pos-10*shift, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-10*shift, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-11*shift, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-11*shift, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            % 5. Quantification
                            text(0.4, text_pos-12*shift, 'Quantification ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');
                            
                            str1 = 'GABA+/Cr: ';
                            str2 = sprintf('%.2f', MRS_struct.out.(vox{kk}).GABA.ConcCr(ii));
                            str3 = 'Glx/Cr: ';
                            str4 = sprintf('%.2f', MRS_struct.out.(vox{kk}).Glx.ConcCr(ii));
                            
                            text(0.4, text_pos-13*shift, str1, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-13*shift, str2, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            text(0.4, text_pos-14*shift, str3, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                            text(0.425, text_pos-14*shift, str4, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                            
                            n = 0;
                            
                    end
                    
                    % 6. FitVer
                    text(0.4, text_pos-15.5*shift+n, 'FitVer: ', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj, 'HorizontalAlignment', 'right');
                    text(0.425, text_pos-15.5*shift+n, MRS_struct.info.version.fit, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 10 - font_size_adj);
                    
                end
                
                % Save output as PDF
                run_count = SavePDF(h, MRS_struct, ii, jj, kk, vox, mfilename, run_count, target{jj});

            catch ME

                fprintf('\n');
                warning('********** An error occurred while fitting %s in dataset: ''%s''. Check data. Skipping to next dataset in batch. **********', target{jj}, MRS_struct.metabfile{1,ii});
                error_report{catch_ind} = strrep(sprintf(['Filename: %s\n\n' getReport(ME,'extended','hyperlinks','off') ...
                    '\n\nVisit https://markmikkelsen.github.io/Gannet-docs/index.html for help.'], MRS_struct.metabfile{1,ii}), '\', '\\');
                catch_ind = catch_ind + 1;

            end % end of load-and-processing loop over datasets

            % Display report if errors occurred
            if ~isempty(error_report{1}) && ii == MRS_struct.p.numScans
                opts = struct('WindowStyle', 'non-modal', 'Interpreter', 'tex');
                for ll = flip(1:size(error_report,2))
                    errordlg(error_report{ll}, sprintf('GannetFit Error Report (%d of %d)', ll, size(error_report,2)), opts);
                end
            end

            if ii == MRS_struct.p.numScans && jj ~= length(target)
                fprintf('\n');
            end

        end

    end

    % Reorder structure
    if isfield(MRS_struct, 'mask')
        if isfield(MRS_struct, 'waterfile')
            structorder = {'info', 'ii', 'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out', 'mask'};
        else
            structorder = {'info', 'ii', 'metabfile', 'p', 'fids', 'spec', 'out', 'mask'};
        end
    else
        if isfield(MRS_struct, 'waterfile')
            structorder = {'info', 'ii', 'metabfile', 'waterfile', 'p', 'fids', 'spec', 'out'};
        else
            structorder = {'info', 'ii', 'metabfile', 'p', 'fids', 'spec', 'out'};
        end
    end
    MRS_struct = orderfields(MRS_struct, structorder);
    
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
        MRS_struct = ExportToCSV(MRS_struct, vox{kk}, 'fit');
    end

end

warning('on','stats:nlinfit:ModelConstantWRTParam');
warning('on','stats:nlinfit:IllConditionedJacobian');
warning('on','stats:nlinfit:IterationLimitExceeded');
warning('on','MATLAB:rankDeficientMatrix');

% Need to close hidden figures to show figures after Gannet is done running
if MRS_struct.p.hide && exist('figTitle','var')
    close(figTitle);
end

end
