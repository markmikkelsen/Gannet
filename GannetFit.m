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
MRS_struct.info.version.fit = '250914';

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
lsqopts = optimset(lsqopts,'MaxIter',800,'TolX',1e-4,'TolFun',1e-4,'Display','off');
nlinopts = statset('nlinfit');
nlinopts = statset(nlinopts,'MaxIter',400,'TolX',1e-6,'TolFun',1e-6,'FunValCheck','off');

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
                    
                    case 'GABA'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   GABA
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            
                            freqbounds = find(freq <= 3.55 & freq >= 2.6);
                            plotbounds = find(freq <= 3.6 & freq >= 2.6);
                            
                            maxinGABA = abs(max(real(DIFF(ii,freqbounds))) - min(real(DIFF(ii,freqbounds))));
                            grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                            LinearInit = grad_points ./ abs(freq(1) - freq(2));
                            constInit = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                            
                            GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
                            GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow
                            
                            lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
                            ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
                            lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
                            ub([1 4 5]) = ub([1 4 5]) / maxinGABA;
                            
                            % Down-weight co-edited Cho signal by including
                            % observation weights in nonlinear regression
                            w = ones(size(DIFF(ii,freqbounds)));
                            residfreq = freq(freqbounds);
                            ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                            weightRange = ChoRange;
                            w(weightRange) = 0.001;
                            
                            % Least-squares model fitting
                            GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, lb, ub, lsqopts);
                            modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
                            [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                            [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                            
                            % Rescale fit parameters and residuals
                            GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
                            resid = resid * maxinGABA;
                            residPlot = residPlot * maxinGABA;
                            
                        else
                            
                            freqbounds = find(freq <= 3.55 & freq >= 2.79);
                            plotbounds = find(freq <= 3.6 & freq >= 2.7);
                            
                            maxinGABA = abs(max(real(DIFF(ii,freqbounds))) - min(real(DIFF(ii,freqbounds))));
                            grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                            LinearInit = grad_points ./ abs(freq(1) - freq(2));
                            constInit = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                            
                            GaussModelInit = [maxinGABA -90 3.026 -LinearInit constInit];
                            GaussModelInit([1 4 5]) = GaussModelInit([1 4 5]) / maxinGABA; % Scale initial conditions to avoid warnings about numerical underflow
                            
                            lb = [0 -200 2.87 -40*maxinGABA -2000*maxinGABA];
                            ub = [4000*maxinGABA -40 3.12 40*maxinGABA 1000*maxinGABA];
                            lb([1 4 5]) = lb([1 4 5]) / maxinGABA;
                            ub([1 4 5]) = ub([1 4 5]) / maxinGABA;
                            
                            % Down-weight Cho subtraction artifact by including
                            % observation weights in nonlinear regression; improves
                            % accuracy of peak fitting (MM: 170701 - thanks to Alex
                            % Craven of University of Bergen for this idea)
                            w = ones(size(DIFF(ii,freqbounds)));
                            residfreq = freq(freqbounds);
                            ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                            weightRange = ChoRange;
                            w(weightRange) = 0.001;
                            
                            % Weighted least-squares model fitting
                            GaussModelInit = lsqcurvefit(@GaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, lb, ub, lsqopts);
                            modelFun_w = @(x,freq) sqrt(w) .* GaussModel(x,freq); % add weights to the model
                            [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGABA, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                            [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGABA, @GaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                            
                            % Rescale fit parameters and residuals
                            GaussModelParam([1 4 5]) = GaussModelParam([1 4 5]) * maxinGABA;
                            resid = resid * maxinGABA;
                            residPlot = residPlot * maxinGABA;
                            
                        end
                        
                        GABAheight = GaussModelParam(1);
                        MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / GABAheight;
                        MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = GaussModelParam(1) ./ sqrt(-GaussModelParam(2)) * sqrt(pi);
                        sigma = sqrt(1/(2*(abs(GaussModelParam(2)))));
                        MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2 * MRS_struct.p.LarmorFreq(ii)) * sigma);
                        MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                        MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;
                        
                        % Calculate SNR of GABA signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GABAheight) / noiseSigma_DIFF;
                        
                    case 'GSH'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   GSH
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        freqbounds = find(freq <= 3.5 & freq >= 2.25);
                        plotbounds = find(freq <= 4.2 & freq >= 1.75);
                        
                        GSHbounds = freq <= 3.3 & freq >= 2.85;
                        Aspbounds = freq <= 2.85 & freq >= 2.25;
                        
                        maxinGSH = max(abs(real(DIFF(ii,GSHbounds))));
                        [maxinAsp, maxInd] = max(abs(real(DIFF(ii,Aspbounds))));
                        
                        offset      = real(DIFF(ii,freqbounds(end)));
                        grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                        LinearInit  = grad_points ./ abs(freq(1) - freq(2));
                        
                        DIFF_Asp = DIFF(ii,Aspbounds);
                        s = sign(real(DIFF_Asp(maxInd)));
                        maxinAsp = s * maxinAsp;
                        
                        if MRS_struct.p.HERMES
                            s = -1;
                        else
                            s = 1;
                        end
                        
                        if MRS_struct.p.TE(ii) < 100
                            
                            GSHgaussModel = @FiveGaussModel;
                            
                            GaussModelInit = [maxinGSH        -300  2.95 ...
                                              s*maxinAsp*0.25 -500  2.73 ...
                                              maxinAsp        -1000 2.61 ...
                                              maxinAsp        -1000 2.55 ...
                                              s*maxinAsp*0.15 -600  2.45 ...
                                              offset -LinearInit -LinearInit];
                            GaussModelInit([1 4 7 10 13 16 17 18]) = GaussModelInit([1 4 7 10 13 16 17 18]) / maxinGSH; % Scale initial conditions to avoid warnings about numerical underflow
                            
                            lb = [-4000*maxinGSH       -1000 2.95-0.02 ...
                                   4000*maxinAsp*0.25  -1000 2.73-0.02 ...
                                   4000*maxinAsp       -1000 2.61-0.02 ...
                                   4000*maxinAsp       -1000 2.55-0.02 ...
                                   4000*maxinAsp*0.15  -1000 2.45-0.02 ...
                                  -2000*abs(offset) 2000*maxinAsp 2000*maxinAsp];
                            ub =  [4000*maxinGSH      -40 2.95+0.02 ...
                                  -4000*maxinAsp*0.25 -40 2.73+0.02 ...
                                  -4000*maxinAsp      -40 2.61+0.02 ...
                                  -4000*maxinAsp      -40 2.55+0.02 ...
                                  -4000*maxinAsp*0.15 -40 2.45+0.02 ...
                                   1000*abs(offset) -1000*maxinAsp -1000*maxinAsp];
                            lb([1 4 7 10 13 16 17 18]) = lb([1 4 7 10 13 16 17 18]) / maxinGSH;
                            ub([1 4 7 10 13 16 17 18]) = ub([1 4 7 10 13 16 17 18]) / maxinGSH;
                            
                        else
                            
                            GSHgaussModel = @SixGaussModel;
                            
                            GaussModelInit = [maxinGSH           -300  2.95 ...
                                              maxinAsp*0.7  -500  2.73 ...
                                              maxinAsp      -1000 2.63 ...
                                              maxinAsp*0.7  -1000 2.58 ...
                                              maxinAsp*0.5  -600  2.46 ...
                                              maxinAsp*0.35 -600  2.37 ...
                                              offset -LinearInit -LinearInit];
                            GaussModelInit([1 4 7 10 13 16 19 20 21]) = GaussModelInit([1 4 7 10 13 16 19 20 21]) / maxinGSH; % Scale initial conditions to avoid warnings about numerical underflow
                            
                            lb = [-4000*maxinGSH           -1000 2.95-0.02 ...
                                  -4000*maxinAsp*0.7  -1000 2.73-0.02 ...
                                  -4000*maxinAsp      -1000 2.63-0.02 ...
                                  -4000*maxinAsp*0.7  -1000 2.58-0.02 ...
                                  -4000*maxinAsp*0.5  -1000 2.46-0.02 ...
                                  -4000*maxinAsp*0.35 -1000 2.37-0.02 ...
                                  -2000*abs(offset) -2000*maxinAsp -2000*maxinAsp];
                            ub =  [4000*maxinGSH           -40 2.95+0.02 ...
                                   4000*maxinAsp*0.7  -40 2.73+0.02 ...
                                   4000*maxinAsp      -40 2.63+0.02 ...
                                   4000*maxinAsp*0.7  -40 2.58+0.02 ...
                                   4000*maxinAsp*0.5  -40 2.46+0.02 ...
                                   4000*maxinAsp*0.35 -40 2.37+0.02 ...
                                   1000*abs(offset) 1000*maxinAsp 1000*maxinAsp];
                            lb([1 4 7 10 13 16 19 20 21]) = lb([1 4 7 10 13 16 19 20 21]) / maxinGSH;
                            ub([1 4 7 10 13 16 19 20 21]) = ub([1 4 7 10 13 16 19 20 21]) / maxinGSH;
                            
                        end
                        
                        if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target, {'EtOH','GABA','GSH'}))
                            w = ones(size(DIFF(ii,freqbounds)));
                            residfreq = freq(freqbounds);
                            ChoRange = residfreq >= 3.13 & residfreq <= 3.3;
                            weightRange = ChoRange;
                            w(weightRange) = 0.001;
                        else
                            w = ones(size(DIFF(ii,freqbounds)));
                        end
                        
                        % Weighted least-squares model fitting
                        GaussModelInit = lsqcurvefit(GSHgaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGSH, lb, ub, lsqopts);
                        modelFun_w = @(x,freq) sqrt(w) .* GSHgaussModel(x,freq); % add weights to the model
                        [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGSH, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                        [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGSH, GSHgaussModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                        
                        % Rescale fit parameters and residuals
                        if MRS_struct.p.TE(ii) < 100
                            GaussModelParam([1 4 7 10 13 16 17 18]) = GaussModelParam([1 4 7 10 13 16 17 18]) * maxinGSH;
                        else
                            GaussModelParam([1 4 7 10 13 16 19 20 21]) = GaussModelParam([1 4 7 10 13 16 19 20 21]) * maxinGSH;
                        end
                        resid = resid * maxinGSH;
                        residPlot = residPlot * maxinGSH;
                        
                        GSHGaussModelParam = GaussModelParam;
                        if MRS_struct.p.TE(ii) < 100 % MM (190613)
                            GSHGaussModelParam(4:3:13) = 0;
                        else
                            GSHGaussModelParam(4:3:16) = 0;
                        end
                        
                        BaselineModelParam = GSHGaussModelParam;
                        BaselineModelParam(1) = 0;
                        
                        MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = real(sum(FiveGaussModel(GSHGaussModelParam, freq(freqbounds)) - FiveGaussModel(BaselineModelParam, freq(freqbounds)))) ...
                            * abs(freq(1) - freq(2));
                        GSHheight = GSHGaussModelParam(1);
                        
                        % Range to determine residuals for GSH
                        residfreq = freq(freqbounds);
                        residGSH  = resid(residfreq <= 3.3 & residfreq >= 2.82);
                        
                        MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(residGSH) / GSHheight;
                        sigma = sqrt(1/(2*(abs(GSHGaussModelParam(2)))));
                        MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) =  abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                        MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                        MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = residGSH;
                        
                        % Calculate SNR of GSH signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).(target{jj}).SNR(ii) = abs(GSHheight) / noiseSigma_DIFF;
                        
                        % MM (200728)
                        %MRS_struct.out.(vox{kk}).(target{jj}).FitError2(ii) = sqrt(mean(residGSH.^2)) / (0.5*noiseSigma_DIFF);
                        
                    case 'Lac'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   Lac
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        freqbounds = find(freq <= 1.8 & freq >= 0.5);
                        plotbounds = find(freq <= 2.12 & freq >= 0);
                        
                        offset   = (mean(real(DIFF(ii,freqbounds(1:10))),2) + mean(real(DIFF(ii,freqbounds((end-9):end))),2)) / 2;
                        slope    = (mean(real(DIFF(ii,freqbounds(1:10))),2) - mean(real(DIFF(ii,freqbounds((end-9):end))),2)) / abs(freq(freqbounds(1)) - freq(freqbounds(end)));
                        maxinLac = max(real(DIFF(ii,freqbounds)));
                        
                        LacModelInit = [maxinLac/2 -1000 1.31 ...
                                        maxinLac/4 -90   1.21 ...
                                        offset slope 0 ...
                                        maxinLac/2];
                        lb = [0 -4000 1.31-0.02 ...
                              0 -100  1.24-0.02 ...
                             -1 -1 -1 ...
                              0];
                        ub = [maxinLac -500 1.31+0.02 ...
                              maxinLac  0   1.24+0.02 ...
                              1 1 1 ...
                              maxinLac];
                        
                        LacModelInit = lsqcurvefit(@LacModel, LacModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub,lsqopts);
                        [LacPlusModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @LacModel, LacModelInit, nlinopts);
                        
                        MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:) = LacPlusModelParam;
                        %MMmodelParam         = LacPlusModelParam;
                        %MMmodelParam([4 7])  = 0;
                        LacModelParam         = LacPlusModelParam;
                        LacModelParam([4 10]) = 0;
                        BHBmodelParam         = LacPlusModelParam;
                        BHBmodelParam(1)      = 0;
                        
                        MRS_struct.out.(vox{kk}).Lac.Area(ii) = sum(LacModel([LacPlusModelParam(1:3) zeros(1,7)], freq(freqbounds))) * abs(freq(1) - freq(2)); % NB: this is Lac+
                        modelHeight = max(LacModel([LacPlusModelParam(1:6) zeros(1,3) LacPlusModelParam(10)], freq(freqbounds)));
                        MRS_struct.out.(vox{kk}).Lac.FitError(ii) = 100 * std(resid) / modelHeight;
                        MRS_struct.out.(vox{kk}).Lac.FWHM(ii) = NaN; % MM (170818): Still need to calculate FWHM
                        MRS_struct.out.(vox{kk}).Lac.Resid(ii,:) = resid;
                        
                        % Calculate SNR of Lac signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).Lac.SNR(ii) = abs(modelHeight) / noiseSigma_DIFF;
                        
                    case 'Glx'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   Glx
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        freqbounds = find(freq <= 4.1 & freq >= 3.45);
                        plotbounds = find(freq <= 4.5 & freq >= 3);
                        
                        maxinGlx    = max(real(DIFF(ii,freqbounds)));
                        grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                        LinearInit  = grad_points ./ abs(freq(1) - freq(2));
                        constInit   = (real(DIFF(ii,freqbounds(end))) + real(DIFF(ii,freqbounds(1))))./2;
                        
                        GaussModelInit = [maxinGlx -90 3.72 maxinGlx -90 3.77 -LinearInit constInit];
                        lb = [0 -200 3.72-0.01 0 -200 3.77-0.01 -40*maxinGlx -2000*maxinGlx];
                        ub = [4000*maxinGlx -40 3.72+0.01 4000*maxinGlx -40 3.77+0.01 40*maxinGlx 1000*maxinGlx];
                        
                        GaussModelInit = lsqcurvefit(@DoubleGaussModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)), lb, ub, lsqopts);
                        [GaussModelParam, resid] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)), @DoubleGaussModel, GaussModelInit, nlinopts);
                        
                        Glxheight = max(GaussModelParam([1,4]));
                        MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / Glxheight;
                        MRS_struct.out.(vox{kk}).(target{jj}).Area(ii)     = (GaussModelParam(1) / sqrt(-GaussModelParam(2)) * sqrt(pi)) + ...
                            (GaussModelParam(4) / sqrt(-GaussModelParam(5)) * sqrt(pi));
                        sigma = (sqrt(1/(2*(abs(GaussModelParam(2)))))) + (sqrt(1/(2*(abs(GaussModelParam(5))))));
                        MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * sigma);
                        MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = GaussModelParam;
                        MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:) = resid;
                        
                        % Calculate SNR of Glx signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight) / noiseSigma_DIFF;
                        
                    case 'GABAGlx'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   GABA+Glx
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        freqbounds = find(freq <= 4.1 & freq >= 2.79);
                        plotbounds = find(freq <= 4.2 & freq >= 2.7);
                        
                        GABAbounds = freq <= 3.2 & freq >= 2.79;
                        Glxbounds  = freq <= 4.1 & freq >= 3.4;
                        
                        maxinGABA   = max(real(DIFF(ii,GABAbounds)));
                        maxinGlx    = max(real(DIFF(ii,Glxbounds)));
                        grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                        LinearInit  = grad_points ./ abs(freq(1) - freq(2));
                        
                        GaussModelInit = [maxinGlx -700 3.71 maxinGlx -700 3.79 maxinGABA -90 3.02 -LinearInit 0 0];
                        GaussModelInit([1 4 7 10]) = GaussModelInit([1 4 7 10]) / maxinGlx; % Scale initial conditions to avoid warnings about numerical underflow
                        
                        lb = [-4000*maxinGlx -1000 3.71-0.02 -4000*maxinGlx -1000 3.79-0.02 -4000*maxinGABA -200 3.02-0.05 -40*maxinGABA -2000*maxinGABA -2000*maxinGABA];
                        ub = [4000*maxinGlx  -40   3.71+0.02  4000*maxinGlx -40   3.79+0.02  4000*maxinGABA -40  3.02+0.05  40*maxinGABA  1000*maxinGABA  1000*maxinGABA];
                        lb([1 4 7 10]) = lb([1 4 7 10]) / maxinGlx;
                        ub([1 4 7 10]) = ub([1 4 7 10]) / maxinGlx;
                        
                        % Down-weight Cho subtraction artifact and (if HERMES)
                        % signals downfield of Glx by including observation weights
                        % in nonlinear regression; improves accuracy of peak
                        % fittings (MM: 170701 - thanks to Alex Craven of
                        % University of Bergen for this idea)
                        w = ones(size(DIFF(ii,freqbounds)));
                        residfreq = freq(freqbounds);
                        ChoRange = residfreq >= 3.16 & residfreq <= 3.285;
                        GlxDownfieldRange = residfreq >= 3.9 & residfreq <= 4.2;
                        if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor, {'Philips','Philips_data','Philips_raw'}))
                            weightRange = ChoRange | GlxDownfieldRange;
                        else
                            weightRange = ChoRange;
                        end
                        w(weightRange) = 0.001;
                        
                        % Weighted least-squares model fitting
                        GaussModelInit = lsqcurvefit(@GABAGlxModel, GaussModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGlx, lb, ub, lsqopts);
                        modelFun_w = @(x,freq) sqrt(w) .* GABAGlxModel(x,freq); % add weights to the model
                        [GaussModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinGlx, modelFun_w, GaussModelInit, nlinopts); % add weights to the data
                        [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinGlx, @GABAGlxModel, GaussModelParam, nlinopts); % re-run for residuals for output figure
                        
                        % Rescale fit parameters and residuals
                        GaussModelParam([1 4 7 10 11 12]) = GaussModelParam([1 4 7 10 11 12]) * maxinGlx;
                        resid     = resid * maxinGlx;
                        residPlot = residPlot * maxinGlx;
                        
                        % Range to determine residuals for GABA and Glx
                        residGABA = resid(residfreq <= 3.55 & residfreq >= 2.79);
                        residGlx  = resid(residfreq <= 4.10 & residfreq >= 3.45);
                        
                        % GABA fitting output
                        MRS_struct.out.(vox{kk}).GABA.Area(ii) = GaussModelParam(7) ./ sqrt(-GaussModelParam(8)) * sqrt(pi);
                        GABAheight = GaussModelParam(7);
                        MRS_struct.out.(vox{kk}).GABA.FitError(ii) = 100 * std(residGABA) / GABAheight;
                        sigma = sqrt(1/(2*(abs(GaussModelParam(8)))));
                        MRS_struct.out.(vox{kk}).GABA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii))*sigma);
                        MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:) = GaussModelParam;
                        MRS_struct.out.(vox{kk}).GABA.Resid(ii,:) = residGABA;
                        
                        % Calculate SNR of GABA signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).GABA.SNR(ii) = abs(GABAheight) / noiseSigma_DIFF;
                        
                        % MM (200728)
                        %MRS_struct.out.(vox{kk}).GABA.FitError2(ii) = sqrt(mean(residGABA.^2)) / (0.5*noiseSigma_DIFF);
                        
                        % Glx fitting output
                        MRS_struct.out.(vox{kk}).Glx.Area(ii) = (GaussModelParam(1) / sqrt(-GaussModelParam(2)) * sqrt(pi)) + ...
                            (GaussModelParam(4) / sqrt(-GaussModelParam(5)) * sqrt(pi));
                        Glxheight = max(GaussModelParam([1,4]));
                        MRS_struct.out.(vox{kk}).Glx.FitError(ii) = 100 * std(residGlx) / Glxheight;
                        sigma = sqrt(1/(2*(abs(GaussModelParam(2))))) + sqrt(1/(2*(abs(GaussModelParam(5)))));
                        MRS_struct.out.(vox{kk}).Glx.FWHM(ii) = abs(2 * MRS_struct.p.LarmorFreq(ii) * sigma);
                        MRS_struct.out.(vox{kk}).Glx.ModelParam(ii,:) = GaussModelParam;
                        MRS_struct.out.(vox{kk}).Glx.Resid(ii,:) = residGlx;
                        
                        % Calculate SNR of Glx signal
                        MRS_struct.out.(vox{kk}).Glx.SNR(ii) = abs(Glxheight) / noiseSigma_DIFF;
                        
                        % MM (200728)
                        %MRS_struct.out.(vox{kk}).Glx.FitError2(ii) = sqrt(mean(residGlx.^2)) / (0.5*noiseSigma_DIFF);
                        
                    case 'EtOH'
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %   EtOH
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        freqbounds = find(freq <= 1.8 & freq >= 0.6);
                        plotbounds = find(freq <= 1.9 & freq >= 0.4);
                        
                        maxinEtOH = max(real(DIFF(ii,freqbounds)));
                        grad_points = (real(DIFF(ii,freqbounds(end))) - real(DIFF(ii,freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1));
                        LinearInit = grad_points ./ abs(freq(1) - freq(2));
                        
                        LorentzModelInit = [maxinEtOH 1.11 1/500 ...
                                            maxinEtOH 1.23 1/500 ...
                                            -LinearInit 0];
                        LorentzModelInit([1 4 7]) = LorentzModelInit([1 4 7]) / maxinEtOH; % Scale initial conditions to avoid warnings about numerical underflow
                        
                        lb = [0             1.11-0.01 1/700 0             1.23-0.01 1/700 -40*maxinEtOH -2e3*maxinEtOH];
                        ub = [100*maxinEtOH 1.11+0.01 1/300 100*maxinEtOH 1.23+0.01 1/300  40*maxinEtOH  1e3*maxinEtOH];
                        
                        lb([1 4 7]) = lb([1 4 7]) / maxinEtOH;
                        ub([1 4 7]) = ub([1 4 7]) / maxinEtOH;
                        
                        % Down-weight co-edited Lac signal by including observation
                        % weights in nonlinear regression
                        w = ones(size(DIFF(ii,freqbounds)));
                        residfreq = freq(freqbounds);
                        LacRange = residfreq >= 1.29 & residfreq <= 1.51;
                        weightRange = LacRange;
                        w(weightRange) = 0.001;
                        
                        % Weighted least-squares model fitting
                        LorentzModelInit = lsqcurvefit(@EtOHModel, LorentzModelInit, freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinEtOH, lb, ub, lsqopts);
                        modelFun_w = @(x,freq) sqrt(w) .* EtOHModel(x,freq); % add weights to the model
                        [LorentzModelParam, resid] = nlinfit(freq(freqbounds), sqrt(w) .* real(DIFF(ii,freqbounds)) / maxinEtOH, modelFun_w, LorentzModelInit, nlinopts); % add weights to the data
                        [~, residPlot] = nlinfit(freq(freqbounds), real(DIFF(ii,freqbounds)) / maxinEtOH, @EtOHModel, LorentzModelParam, nlinopts); % re-run for residuals for output figure
                        
                        % Rescale fit parameters and residuals
                        LorentzModelParam([1 4 7 8]) = LorentzModelParam([1 4 7 8]) * maxinEtOH;
                        resid = resid * maxinEtOH;
                        residPlot = residPlot * maxinEtOH;
                        
                        EtOHheight = max(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqbounds)));
                        MRS_struct.out.(vox{kk}).(target{jj}).FitError(ii) = 100 * std(resid) / EtOHheight;
                        MRS_struct.out.(vox{kk}).(target{jj}).Area(ii) = sum(EtOHModel([LorentzModelParam(1:end-2) 0 0],freq(freqbounds))) * abs(freq(1) - freq(2));
                        
                        MRS_struct.out.(vox{kk}).(target{jj}).FWHM(ii)         = (LorentzModelParam(3) + LorentzModelParam(6)) * MRS_struct.p.LarmorFreq(ii);
                        MRS_struct.out.(vox{kk}).(target{jj}).ModelParam(ii,:) = LorentzModelParam;
                        MRS_struct.out.(vox{kk}).(target{jj}).Resid(ii,:)      = resid;
                        
                        % Calculate SNR of EtOH signal
                        noiseSigma_DIFF = CalcNoise(freq, DIFF(ii,:));
                        MRS_struct.out.(vox{kk}).EtOH.SNR(ii) = abs(EtOHheight) / noiseSigma_DIFF;
                        
                    otherwise
                        
                        error('Metabolite ''%s'' not recognized.', target{jj});
                        
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   1a. Initialize the output figure
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
                fig_w = 1000;
                fig_h = 707;
                set(h,'Position',[(scr_sz(3)-fig_w)/2, (scr_sz(4)-fig_h)/2, fig_w, fig_h]);
                set(h,'Color',[1 1 1]);
                figTitle = 'GannetFit Output';
                set(h,'Name',figTitle,'Tag',figTitle,'NumberTitle','off');
                
                % Spectra plot
                subplot(2,2,1);
                metabmin = min(real(DIFF(ii,plotbounds)));
                metabmax = max(real(DIFF(ii,plotbounds)));
                resmax = max(resid);
                resid = resid + metabmin - resmax;
                
                switch target{jj}
                    case 'GABA'
                        residPlot = residPlot + metabmin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(weightRange) = NaN;
                        hold on;
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                            freq(freqbounds), GaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), residPlot2, 'k');
                        plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                        hold off;
                        set(gca, 'XLim', [2.6 3.6], 'FontSize', 10 - font_size_adj);
                        
                    case 'GSH'
                        residPlot = residPlot + metabmin - max(residPlot);
                        residPlot2 = residPlot;
                        if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            residPlot2(weightRange) = NaN;
                            hold on;
                            plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b' , ...
                                freq(freqbounds), GSHgaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                                freq(freqbounds), residPlot2, 'k');
                            plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                            hold off;
                        else
                            plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b' , ...
                                freq(freqbounds), GSHgaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                                freq(freqbounds), residPlot2, 'k');
                        end
                        set(gca, 'XLim', [1.8 4.2], 'FontSize', 10 - font_size_adj);
                        
                    case 'Lac'
                        hold on;
                        p1 = plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'k');
                        p2 = plot(freq(freqbounds), LacModel(BHBmodelParam,freq(freqbounds)), 'Color', [31 120 180]/255, 'LineWidth', 1);
                        p3 = plot(freq(freqbounds), LacModel(LacModelParam,freq(freqbounds)), 'Color', [51 160 44]/255, 'LineWidth', 1);
                        p4 = plot(freq(freqbounds), LacModel(LacPlusModelParam,freq(freqbounds)), 'r', 'LineWidth', 1);
                        p5 = plot(freq(freqbounds), resid, 'k');
                        hold off;
                        set(gca, 'XLim', [0.7 1.9], 'XTick', 0:0.25:10, 'FontSize', 10 - font_size_adj);
                        
                    case 'Glx'
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                            freq(freqbounds), DoubleGaussModel(GaussModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), resid, 'k');
                        set(gca, 'XLim', [3.4 4.2], 'FontSize', 10 - font_size_adj);
                        
                    case 'GABAGlx'
                        residPlot  = residPlot + metabmin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(weightRange) = NaN;
                        hold on;
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                            freq(freqbounds), GABAGlxModel(GaussModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), residPlot2, 'k');
                        % Plot weighted portion of residuals in different color
                        if MRS_struct.p.HERMES && any(strcmp(MRS_struct.p.vendor,{'Philips','Philips_data','Philips_raw'}))
                            plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                            plot(freq(freqbounds(GlxDownfieldRange)), residPlot(GlxDownfieldRange), 'Color', [255 160 64]/255);
                        else
                            plot(freq(freqbounds(ChoRange)), residPlot(ChoRange), 'Color', [255 160 64]/255);
                        end
                        hold off;
                        set(gca, 'XLim', [2.7 4.2], 'FontSize', 10 - font_size_adj);
                        
                    case 'EtOH'
                        residPlot = residPlot + metabmin - max(residPlot);
                        residPlot2 = residPlot;
                        residPlot2(weightRange) = NaN;
                        hold on;
                        plot(freq(plotbounds), real(DIFF(ii,plotbounds)), 'b', ...
                            freq(freqbounds), EtOHModel(LorentzModelParam,freq(freqbounds)), 'r', ...
                            freq(freqbounds), residPlot2, 'k');
                        plot(freq(freqbounds(LacRange)), residPlot(LacRange), 'Color', [255 160 64]/255);
                        hold off;
                        set(gca, 'XLim', [0.4 1.9], 'FontSize', 10 - font_size_adj);
                end
                
                % From here on is cosmetic - adding labels etc.
                switch target{jj}
                    case 'GABA'
                        text(3, metabmax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelbounds = freq <= 2.4 & freq >= 2;
                        tailtop = max(real(DIFF(ii,labelbounds)));
                        tailbottom = min(real(DIFF(ii,labelbounds)));
                        text(2.8, min(resid), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(2.8, tailtop + metabmax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.8, tailbottom - metabmax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        if length(MRS_struct.p.target) ~= 3
                            text(3.2, min(residPlot) - 0.5*abs(max(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                                'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        end
                        
                    case 'GSH'
                        text(2.95, maxinGSH+maxinGSH/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelbounds = freq <= 2.4 & freq >= 1.75;
                        tailtop = max(real(DIFF(ii,labelbounds)));
                        tailbottom = min(real(DIFF(ii,labelbounds)));
                        text(2.25, min(resid), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(2.25, tailtop + metabmax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.45, tailbottom - 20*metabmax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        if length(MRS_struct.p.target) == 3 && all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            text(3.3, max(residPlot) + 0.5*abs(max(residPlot) - min(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                                'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
                        end
                        
                    case 'Lac'
%                         text(1.27, metabmax/3.5, 'Lac+', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
%                         labelbounds = freq <= 0.8 & freq >= 0.42;
%                         tailtop = max(real(DIFF(ii,labelbounds)));
%                         tailbottom = min(real(DIFF(ii,labelbounds)));
%                         text(0.44, mean(real(DIFF(ii,labelbounds))) + 4*std(real(DIFF(ii,labelbounds))), 'data', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
%                         text(0.44, mean(resid) - 4*std(resid), 'residual', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
%                         text(0.85, tailbottom, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        legend([p1, p2, p3, p4, p5], {'data','BHB+','Lac+','model','residual'}, 'box', 'off', 'Location', 'northeast', 'FontSize', 9 - font_size_adj);
                        
                    case 'Glx'
                        text(3.8, metabmax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelbounds = freq <= 3.6 & freq >= 3.4;
                        tailtop = max(real(DIFF(ii,labelbounds)));
                        tailbottom = min(real(DIFF(ii,labelbounds)));
                        text(3.5, min(resid),'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(3.5, tailtop + metabmax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(3.5, tailbottom - metabmax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        
                    case 'GABAGlx'
                        text(3, metabmax/4, 'GABA', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        text(3.755, metabmax/4, 'Glx', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        text(2.8, min(resid), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        labelbounds = freq <= 2.8 & freq >= 2.7;
                        tailtop = max(real(DIFF(ii,labelbounds)));
                        tailbottom = min(real(DIFF(ii,labelbounds)));
                        text(2.8, tailtop + metabmax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(2.8, tailbottom - metabmax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                        text(3.2, min(residPlot) - 0.5*abs(max(residPlot)), 'down-weighted', 'Color', [255 160 64]/255, ...
                            'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        
                    case 'EtOH'
                        text(1.45, metabmax/4, target{jj}, 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                        labelbounds = freq <= 0.9 & freq >= 0.5;
                        tailtop = max(real(DIFF(ii,labelbounds)));
                        tailbottom = min(real(DIFF(ii,labelbounds)));
                        text(0.6, min(resid),'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                        text(0.8, tailtop + metabmax/20, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                        text(0.8, tailbottom - metabmax/20, 'model', 'Color', [1 0 0], 'FontSize', 10 - font_size_adj);
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
                    
                    % Estimate height and baseline from data
                    [maxinWater, watermaxindex] = max(real(WaterData(ii,:)),[],2);
                    waterbase = mean(real(WaterData(ii,freq <= 4 & freq >= 3.8)));
                    
                    % Philips data do not phase well based on first point, so do a preliminary
                    % fit, then adjust phase of WaterData accordingly
                    % Do this only if the eddy current correction is not performed on water
                    if strcmp(MRS_struct.p.vendor,'Philips') && ~MRS_struct.p.water_ECC
                        % Run preliminary fit of data
                        LGModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50];
                        freqbounds = freq <= 5.6 & freq >= 3.8;
                        
                        % Do the water fit (Lorentz-Gauss)
                        LGModelParam = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModel, LGModelInit, nlinopts);
                        
                        Eerror = zeros([120 1]);
                        for ll = 1:120
                            Data = WaterData(ii,freqbounds)*exp(1i*pi/180*ll*3);
                            Model = LorentzGaussModel(LGModelParam,freq(freqbounds));
                            Eerror(ll) = sum((real(Data)-Model).^2);
                        end
                        [~,index] = min(Eerror);
                        WaterData(ii,:) = WaterData(ii,:) * exp(1i*pi/180*index*3);
                    end
                    
                    LGPModelInit = [maxinWater 20 freq(watermaxindex) 0 waterbase -50 0];
                    lb = [0.01*maxinWater 1   4.6 0     0 -50 -pi];
                    ub = [40*maxinWater   100 4.8 1e-6  1  0   pi];
                    
                    freqbounds = freq <= 5.6 & freq >= 3.8;
                    
                    % Least-squares model fitting
                    LGPModelInit = lsqcurvefit(@LorentzGaussModelP, LGPModelInit, freq(freqbounds), real(WaterData(ii,freqbounds)), lb, ub, lsqopts);
                    [LGPModelParam, residw] = nlinfit(freq(freqbounds), real(WaterData(ii,freqbounds)), @LorentzGaussModelP, LGPModelInit, nlinopts);
                    
                    WaterArea = sum(real(LorentzGaussModel(LGPModelParam(1:end-1), freq(freqbounds))) - BaselineModel(LGPModelParam(3:5), freq(freqbounds)),2);
                    MRS_struct.out.(vox{kk}).water.Area(ii) = WaterArea * abs(freq(1) - freq(2));
                    waterheight = LGPModelParam(1);
                    MRS_struct.out.(vox{kk}).water.FitError(ii) = 100 * std(residw) / waterheight;
                    
                    LG = real(LorentzGaussModel(LGPModelParam(1:end-1), freq(freqbounds))) - BaselineModel(LGPModelParam(3:5), freq(freqbounds));
                    LG = LG/max(LG);
                    ind = find(LG >= 0.5);
                    f = freq(freqbounds);
                    w = abs(f(ind(1)) - f(ind(end)));
                    MRS_struct.out.(vox{kk}).water.FWHM(ii) = w * MRS_struct.p.LarmorFreq(ii);
                    MRS_struct.out.(vox{kk}).water.ModelParam(ii,:) = LGPModelParam;
                    MRS_struct.out.(vox{kk}).water.Resid(ii,:) = residw;
                    
                    % Calculate SNR of water signal
                    noiseSigma_Water = CalcNoise(freq, WaterData(ii,:));
                    MRS_struct.out.(vox{kk}).water.SNR(ii) = abs(waterheight) / noiseSigma_Water;
                    
                    % Water spectrum plot
                    hb = subplot(2,2,3);
                    watmin = min(real(WaterData(ii,:)));
                    watmax = max(real(WaterData(ii,:)));
                    residw = residw + watmin - max(residw);
                    plot(freq(freqbounds), real(WaterData(ii,freqbounds)), 'b', ...
                        freq(freqbounds), real(LorentzGaussModelP(LGPModelParam,freq(freqbounds))), 'r', ...
                        freq(freqbounds), residw, 'k');
                    set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XTick', 4.2:0.2:5.2, 'FontSize', 10 - font_size_adj);
                    xlim([4.2 5.2]);
                    set(get(gca,'YAxis'),'Visible','off');
                    % Add on some labels
                    text(4.8, watmax/2, 'Water', 'HorizontalAlignment', 'right', 'FontSize', 10 - font_size_adj);
                    labelfreq = freq(freqbounds);
                    rlabelbounds = labelfreq <= 4.4 & labelfreq >= 4.25;
                    axis_bottom = axis;
                    text(4.4, max(min(residw(rlabelbounds)) - 0.05 * watmax, axis_bottom(3)), 'residual', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj);
                    xlabel('ppm', 'FontSize', 11 - font_size_adj);
                    title('Reference signals', 'FontSize', 11 - font_size_adj);
                    
                    % Root sum square fit error and concentration in institutional units
                    switch target{jj}
                        case 'GABA'
                            MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);
                            
                        case 'Glx'
                            MRS_struct.out.(vox{kk}).Glx.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);
                            
                        case 'GABAGlx'
                            MRS_struct.out.(vox{kk}).GABA.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GABA.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct.out.(vox{kk}).Glx.FitError_W(ii)  = sqrt(MRS_struct.out.(vox{kk}).Glx.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, 'GABA', ii);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, 'Glx', ii);
                            
                        case 'GSH'
                            MRS_struct.out.(vox{kk}).GSH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).GSH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                            
                        case 'Lac'
                            MRS_struct.out.(vox{kk}).Lac.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).Lac.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                            
                        case 'EtOH'
                            MRS_struct.out.(vox{kk}).EtOH.FitError_W(ii) = sqrt(MRS_struct.out.(vox{kk}).EtOH.FitError(ii).^2 + MRS_struct.out.(vox{kk}).water.FitError(ii).^2);
                            MRS_struct = CalcIU(MRS_struct, vox{kk}, (target{jj}), ii);
                    end
                    
                    % Generate scaled spectra (for plotting)
                    MRS_struct.spec.(vox{kk}).(target{jj}).off_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    if MRS_struct.p.HERMES
                        MRS_struct.spec.(vox{kk}).(target{jj}).off_off_scaled(ii,:) = ...
                            MRS_struct.spec.(vox{kk}).(target{jj}).off_off(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    end
                    MRS_struct.spec.(vox{kk}).(target{jj}).on_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).on(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    MRS_struct.spec.(vox{kk}).(target{jj}).diff_unfilt_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).diff_unfilt(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    MRS_struct.spec.(vox{kk}).(target{jj}).diff_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).diff(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    MRS_struct.spec.(vox{kk}).(target{jj}).sum_scaled(ii,:) = ...
                        MRS_struct.spec.(vox{kk}).(target{jj}).sum(ii,:) .* (1/MRS_struct.out.(vox{kk}).water.ModelParam(ii,1));
                    
                    % Reorder structure fields
                    MRS_struct.out.(vox{kk}).water = orderfields(MRS_struct.out.(vox{kk}).water, {'ModelParam', 'Resid', 'Area', 'FWHM', 'SNR', 'FitError'});
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %   2a.  Residual Water Fit
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if MRS_struct.p.fit_resid_water
                        
                        residWater = OFF(ii,:);
                        freqWater  = freq <= 4.68 + 0.4 & freq >= 4.68 - 0.4;
                        residWater = residWater(freqWater);
                        freqWater  = freq(freqWater);
                        
                        [maxResidWater, maxInd] = max(abs(real(residWater)));
                        s = sign(real(residWater(maxInd)));
                        maxResidWater = s * maxResidWater;
                        offset = real(residWater(1));
                        
                        LGPModelInit = [maxResidWater 25 freqWater(maxInd) 0 offset 0.001 0];
                        LGPModelInit([1 5]) = LGPModelInit([1 5]) / maxResidWater; % Scale initial conditions to avoid warnings about numerical underflow
                        
                        %lb = [maxResidWater-abs(2*maxResidWater) 1 freqWaterOFF(maxInd)-0.2 0 0 -200 -pi];
                        %ub = [maxResidWater+abs(2*maxResidWater) 100 freqWaterOFF(maxInd)+0.2 0.000001 1 0 pi];
                        %lb([1 4 5]) = lb([1 4 5]) / maxResidWater;
                        %ub([1 4 5]) = ub([1 4 5]) / maxResidWater;
                        
                        % Least-squares model fitting
                        [MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,:), residRW] = nlinfit(freqWater, real(residWater) / maxResidWater, @LorentzGaussModelP, LGPModelInit, nlinopts);
                        
                        % Rescale fit parameters and residuals
                        MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) = MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,[1 4 5]) * maxResidWater;
                        residRW = residRW * maxResidWater;
                        
                        MRS_struct.out.(vox{kk}).resid_water.FitError(ii) = 100*std(residRW) / MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1);
                        
                        MRS_struct.out.(vox{kk}).resid_water.SuppressionFactor(ii) = ...
                            (MRS_struct.out.(vox{kk}).water.ModelParam(ii,1) - abs(MRS_struct.out.(vox{kk}).resid_water.ModelParam(ii,1))) ...
                            / MRS_struct.out.(vox{kk}).water.ModelParam(ii,1);
                        
                    end
                    
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   3.  Cr Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Cr_OFF = OFF(ii,:);
                freqboundsChoCr = freq <= 3.6 & freq >= 2.6;
                
                ChoCrMeanSpec   = Cr_OFF(freqboundsChoCr).';
                Baseline_offset = real(ChoCrMeanSpec(1) + ChoCrMeanSpec(end)) / 2;
                Width_estimate  = 0.05;
                Area_estimate   = (max(real(ChoCrMeanSpec)) - min(real(ChoCrMeanSpec))) * Width_estimate * 4;
                ChoCr_initx     = [Area_estimate Width_estimate 3.02 0 Baseline_offset 0 1] ...
                    .* [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
                [ChoCrModelParam, ~, residChoCr] = FitChoCr(freq(freqboundsChoCr), ChoCrMeanSpec, ChoCr_initx, MRS_struct.p.LarmorFreq(ii));
                MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:) = ChoCrModelParam ./ [1 2*MRS_struct.p.LarmorFreq(ii) MRS_struct.p.LarmorFreq(ii) 180/pi 1 1 1];
                
                % Initialise fitting pars
                freqboundsCr = freq <= 3.12 & freq >= 2.72;
                LorentzModelInit = [max(real(Cr_OFF(freqboundsCr))) 0.05 3.0 0 0 0];
                
                % Least-squares model fitting
                LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqboundsCr), real(Cr_OFF(freqboundsCr)), [], [], lsqopts);
                [LorentzModelParam, residCr] = nlinfit(freq(freqboundsCr), real(Cr_OFF(freqboundsCr)), @LorentzModel, LorentzModelInit, nlinopts);
                
                MRS_struct.out.(vox{kk}).Cr.ModelParam(ii,:) = LorentzModelParam;
                CrHeight = LorentzModelParam(1) / (2*pi*LorentzModelParam(2));
                MRS_struct.out.(vox{kk}).Cr.FitError(ii)     = 100 * std(residCr) / CrHeight;
                MRS_struct.out.(vox{kk}).Cr.Resid(ii,:)      = residCr;
                
                MRS_struct.out.(vox{kk}).Cho.ModelParam(ii,:) = MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:);
                ChoHeight = (MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1) * MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,7)) / (2*pi*MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2));
                MRS_struct.out.(vox{kk}).Cho.FitError(ii)     = 100 * std(residChoCr) / ChoHeight;
                MRS_struct.out.(vox{kk}).Cho.Resid(ii,:)      = residChoCr;
                
                MRS_struct.out.(vox{kk}).Cr.Area(ii)  = sum(real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqboundsChoCr)) - ...
                    TwoLorentzModel([0 MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,2:end-1) 0], freq(freqboundsChoCr)))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).Cho.Area(ii) = sum(real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:), freq(freqboundsChoCr)) - ...
                    TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:end-1) 0], freq(freqboundsChoCr)))) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).Cr.FWHM(ii)  = ChoCrModelParam(2);
                MRS_struct.out.(vox{kk}).Cho.FWHM(ii) = ChoCrModelParam(2);
                
                % Calculate SNR of Cr signal
                noiseSigma_OFF = CalcNoise(freq, OFF(ii,:));
                MRS_struct.out.(vox{kk}).Cr.SNR(ii)  = abs(CrHeight) / noiseSigma_OFF;
                MRS_struct.out.(vox{kk}).Cho.SNR(ii) = abs(ChoHeight) / noiseSigma_OFF;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   4.  NAA Fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                NAA_OFF = OFF(ii,:);
                freqbounds = find(freq <= 2.25 & freq >= 1.75);
                
                maxinNAA    = max(real(NAA_OFF(freqbounds)));
                grad_points = (real(NAA_OFF(freqbounds(end))) - real(NAA_OFF(freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1)); %in points
                LinearInit  = grad_points ./ abs(freq(1) - freq(2));
                constInit   = (real(NAA_OFF(freqbounds(end))) + real(NAA_OFF(freqbounds(1)))) ./ 2;
                
                LorentzModelInit = [maxinNAA 0.05 2.01 0 -LinearInit constInit];
                lb = [0 0.01 1.97 0 -40*maxinNAA -2000*maxinNAA];
                ub = [4000*maxinNAA 0.1 2.05 0.5 40*maxinNAA 1000*maxinNAA];
                
                % Least-squares model fitting
                LorentzModelInit = lsqcurvefit(@LorentzModel, LorentzModelInit, freq(freqbounds), real(NAA_OFF(freqbounds)), lb, ub, lsqopts);
                [LorentzModelParam, resid] = nlinfit(freq(freqbounds), real(NAA_OFF(freqbounds)), @LorentzModel, LorentzModelInit, nlinopts);
                
                NAAheight = LorentzModelParam(1) / (2*pi*LorentzModelParam(2));
                MRS_struct.out.(vox{kk}).NAA.FitError(ii) = 100 * std(resid) / NAAheight;
                NAAModelParam = LorentzModelParam;
                NAAModelParam(4) = 0;
                MRS_struct.out.(vox{kk}).NAA.Area(ii) = sum(LorentzModel(NAAModelParam,freq(freqbounds)) - BaselineModel(NAAModelParam([3 6 5]),freq(freqbounds)), 2) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).NAA.FWHM(ii) = abs((2*MRS_struct.p.LarmorFreq(ii)) * NAAModelParam(2));
                MRS_struct.out.(vox{kk}).NAA.ModelParam(ii,:) = LorentzModelParam;
                MRS_struct.out.(vox{kk}).NAA.Resid(ii,:) = resid;

                % Calculate SNR of NAA signal
                MRS_struct.out.(vox{kk}).NAA.SNR(ii) = abs(NAAheight) / noiseSigma_OFF;


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   5.  Glu Fit                       MM (240328): This still needs testing to see if it's reliable
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                Glu_SUM = SUM(ii,:);
                freqbounds = find(freq <= 2.43 & freq >= 2.25);

                maxinGlu    = max(real(Glu_SUM(freqbounds)));
                grad_points = (real(Glu_SUM(freqbounds(end))) - real(Glu_SUM(freqbounds(1)))) ./ abs(freqbounds(end) - freqbounds(1)); %in points
                LinearInit  = grad_points ./ abs(freq(1) - freq(2));
                constInit   = (real(Glu_SUM(freqbounds(end))) + real(Glu_SUM(freqbounds(1)))) ./ 2;

                LorentzModelInit = [maxinGlu -0.2*maxinGlu -0.2*maxinGlu ...
                                    30 ...
                                    2.34 ...
                                    0.05 ...
                                    0 ...
                                    -LinearInit -LinearInit constInit];

                lb = [-4e3*maxinGlu -800*maxinGlu -800*maxinGlu ...
                      0 ...
                      2.34-0.05 ...
                      0 ...
                      -180 ...
                      -40*maxinGlu -40*maxinGlu -2000*maxinGlu];

                ub = [4e3*maxinGlu 800*maxinGlu 800*maxinGlu ...
                      100 ...
                      2.34+0.05 ...
                      0.1 ...
                      180 ...
                      40*maxinGlu 40*maxinGlu 1000*maxinGlu];

                % Least-squares model fitting
                LorentzModelInit = lsqcurvefit(@ThreeLorentzModel, LorentzModelInit, freq(freqbounds), real(Glu_SUM(freqbounds)), lb, ub, lsqopts);
                [LorentzModelParam, resid] = nlinfit(freq(freqbounds), real(Glu_SUM(freqbounds)), @ThreeLorentzModel, LorentzModelInit, nlinopts);

                GluModelParam       = LorentzModelParam;
                GluModelParam(7:10) = 0;

                GluHeight = LorentzModelParam(1) * LorentzModelParam(4);
                MRS_struct.out.(vox{kk}).Glu.FitError(ii)     = 100 * std(resid) / GluHeight;
                MRS_struct.out.(vox{kk}).Glu.Area(ii)         = sum(abs(ThreeLorentzModel(GluModelParam, freq(freqbounds))),2) * abs(freq(1) - freq(2));
                MRS_struct.out.(vox{kk}).Glu.FWHM(ii)         = LorentzModelParam(4) / pi;
                MRS_struct.out.(vox{kk}).Glu.ModelParam(ii,:) = LorentzModelParam;
                MRS_struct.out.(vox{kk}).Glu.Resid(ii,:)      = resid;

                % Calculate SNR of Glu signal
                noiseSigma_Glu_SUM = CalcNoise(freq, SUM(ii,:));
                MRS_struct.out.(vox{kk}).Glu.SNR(ii) = abs(GluHeight) / noiseSigma_Glu_SUM;

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

                Crmin = min(real(Cr_OFF(freqboundsCr)));
                Crmax = max(real(Cr_OFF(freqboundsCr)));
                resmaxCr = max(residCr);
                residCr = residCr + Crmin - resmaxCr;
                if strcmp(MRS_struct.p.reference,'H2O')
                    hb_pos = get(hb, 'Position');
                    hd = axes(hb.Parent, 'Position', [1.3*hb_pos(1)+hb_pos(3)-0.175, 1.001*hb_pos(2)+hb_pos(4)-0.175 0.1125 0.1575]);
                    plot(freq, real(Cr_OFF), 'b', ...
                        freq(freqboundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqboundsChoCr))), 'r', ...
                        freq(freqboundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqboundsChoCr))), 'r', ...
                        freq(freqboundsCr), residCr, 'k');
                    text(2.94, Crmax*0.75, 'Creatine', 'HorizontalAlignment', 'left', 'FontSize', 10 - font_size_adj2);
                    set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XLim', [2.6 3.6], 'XTick', 2.6:0.2:3.6, 'FontSize', 6 - font_size_adj2);
                    xlabel('ppm', 'FontSize', 11 - font_size_adj);
                    set(get(gca,'YAxis'), 'Visible', 'off');
                    hd_pos = get(hd, 'Position');
                    set(hd, 'Position', [hb_pos(1)+hb_pos(3)-hd_pos(3) hd_pos(2:4)]);
                else
                    hb = subplot(2,2,3);
                    plot(freq, real(Cr_OFF), 'b', ...
                        freq(freqboundsChoCr), real(TwoLorentzModel(MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,:),freq(freqboundsChoCr))), 'r', ...
                        freq(freqboundsChoCr), real(TwoLorentzModel([MRS_struct.out.(vox{kk}).ChoCr.ModelParam(ii,1:(end-1)) 0],freq(freqboundsChoCr))), 'r', ...
                        freq(freqboundsCr), residCr, 'k');
                    set(gca,'XDir', 'reverse', 'TickDir', 'out', 'Box', 'off', 'XTick', 2.6:0.2:3.6, 'FontSize', 10 - font_size_adj);
                    set(get(gca,'YAxis'),'Visible','off');
                    xlim([2.6 3.6]);
                    Crlabelbounds = freq(freqboundsCr) <= 3.12 & freq(freqboundsCr) >= 2.72;
                    text(3, max(residCr(Crlabelbounds)) + 0.05*Crmax, 'residual', 'HorizontalAlignment', 'center', 'FontSize', 10 - font_size_adj);
                    text(2.7, 0.1*Crmax, 'data', 'Color', [0 0 1], 'FontSize', 10 - font_size_adj);
                    text(2.7, 0.01*Crmax, 'model','Color', [1 0 0], 'FontSize', 10 - font_size_adj);
                    text(2.94, Crmax*0.75, 'Creatine', 'FontSize', 10 - font_size_adj);
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
