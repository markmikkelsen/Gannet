function PlotPrePostAlign2(MRS_struct, vox, ii)
% Plot pre-/post-alignment spectra
% Updates by MGSaleh 2016, MM 2017-2020

for kk = 1:length(vox)
    
    if MRS_struct.p.HERMES
        
        spectraToPlot = zeros(length(MRS_struct.p.target), length(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff));
        for jj = 1:length(MRS_struct.p.target)
            spectraToPlot(jj,:) = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:);
        end
        
        model      = cell(1,2);
        freqBounds = cell(1,2);
        for jj = 1:length(MRS_struct.p.target)
            switch MRS_struct.p.target{jj}
                case 'GABA'
                    freqBounds{jj} = MRS_struct.spec.freq <= 3.55 & MRS_struct.spec.freq >= 2.79;
                    model{jj}      = GaussModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds{jj}));
                case 'GSH'
                    freqBounds{jj} = MRS_struct.spec.freq <= 3.3 & MRS_struct.spec.freq >= 2.35;
                    model{jj}      = FiveGaussModel(MRS_struct.out.(vox{kk}).GSH.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds{jj}));
                case 'GABAGlx'
                    freqBounds{jj} = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
                    model{jj}      = GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds{jj}));
                case 'Lac'
                    freqBounds{jj} = MRS_struct.spec.freq <= 1.8 & MRS_struct.spec.freq >= 0.5;
                    model{jj}      = LacModel(MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds{jj}));
            end
        end
        
        % Shift baselines to zero
        baseRange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
        % Some bandwidth-limited acquisitions may not record anything below
        % 0 ppm, in this case get the baseline from the other side of water
        if sum(baseRange) == 0
            baseRange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
        end
        peakRange = zeros(length(MRS_struct.p.target), length(MRS_struct.spec.freq));
        for jj = 1:length(MRS_struct.p.target)
            switch MRS_struct.p.target{jj}
                case {'GABA','Glx','GABAGlx'}
                    peakRange(jj,:) = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.2;
                case 'GSH'
                    peakRange(jj,:) = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
                case {'Lac','EtOH'}
                    peakRange(jj,:) = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
            end
        end
        
        specBaseline  = mean(real(spectraToPlot(:,baseRange)),2);
        spectraToPlot = spectraToPlot - repmat(specBaseline, [1 length(spectraToPlot)]);
        for jj = 1:length(MRS_struct.p.target)
            model{jj} = model{jj} - specBaseline(jj);
        end
        
        % Stack spectra
        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'}))
            maxGABAGlx         = abs(max(max(real(spectraToPlot(1,logical(peakRange(1,:)))),[],2)));
            minGSH             = abs(min(min(real(spectraToPlot(2,logical(peakRange(2,:)))),[],2)));
            shift              = max([maxGABAGlx minGSH]) + 0.5*min([maxGABAGlx minGSH]);
            spectraToPlot(2,:) = spectraToPlot(2,:) + shift;
            model{2}           = model{2} + shift;
            
            yAxisMax = max(real(spectraToPlot(2,logical(peakRange(2,:)))));
            yAxisMin = min(real(spectraToPlot(1,logical(peakRange(1,:)))));
            yRange   = abs(yAxisMax - yAxisMin);
            yAxisMax = yAxisMax + 0.1*yRange;
            yAxisMin = yAxisMin - 0.1*yRange;
        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            shift              = abs(min(real(spectraToPlot(3,logical(peakRange(3,:))))));
            spectraToPlot(3,:) = spectraToPlot(3,:) + 1.5*shift;
            shift              = abs(max(real(spectraToPlot(1,logical(peakRange(1,:))))));
            spectraToPlot(1,:) = spectraToPlot(1,:) - 1.5*shift;
            
            yAxisMax = abs(max(real(spectraToPlot(3,logical(peakRange(3,:))))));
            yAxisMax = yAxisMax + 0.2*yAxisMax;
            yAxisMin = max(real(spectraToPlot(1,logical(peakRange(1,:)))));
            yAxisMin = yAxisMin - 4*abs(yAxisMin);
        end
        
        for jj = 1:length(MRS_struct.p.target)
            hold on;
            plot(MRS_struct.spec.freq, real(spectraToPlot(jj,:)), 'Color', 'k');
            plot(MRS_struct.spec.freq(freqBounds{jj}), model{jj}, 'Color', 'r');
            hold off;
        end
        
    else
        
        spectraToPlot = MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:);
        
        switch MRS_struct.p.target{1}
            case 'GABA'
                freqBounds = MRS_struct.spec.freq <= 3.55 & MRS_struct.spec.freq >= 2.79;
                model      = GaussModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds));
            case 'GSH'
                freqBounds = MRS_struct.spec.freq <= 3.3 & MRS_struct.spec.freq >= 2.35;
                model      = FiveGaussModel(MRS_struct.out.(vox{kk}).GSH.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds));
            case 'GABAGlx'
                freqBounds = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.79;
                model      = GABAGlxModel(MRS_struct.out.(vox{kk}).GABA.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds));
            case 'Lac'
                freqBounds = MRS_struct.spec.freq <= 1.8 & MRS_struct.spec.freq >= 0.5;
                model      = LacModel(MRS_struct.out.(vox{kk}).Lac.ModelParam(ii,:),MRS_struct.spec.freq(freqBounds));
        end
        
        % Shift baselines to zero
        baseRange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
        % Some bandwidth-limited acquisitions may not record anything below
        % 0 ppm, in this case get the baseline from the other side of water
        if sum(baseRange) == 0
            baseRange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
        end
        switch MRS_struct.p.target{1}
            case {'GABA','Glx','GABAGlx'}
                peakRange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.26;
            case 'GSH'
                peakRange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
            case {'Lac','EtOH'}
                peakRange = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
        end
        
        specBaseline  = mean(real(spectraToPlot(baseRange)));
        spectraToPlot = spectraToPlot - specBaseline;
        model         = model - specBaseline;
        
        hold on;
        plot(MRS_struct.spec.freq, real(spectraToPlot), 'Color', 'k');
        plot(MRS_struct.spec.freq(freqBounds), model, 'Color', 'r');
        hold on;
        
        yAxisMax = max(real(spectraToPlot(peakRange)));
        yAxisMin = min(real(spectraToPlot(peakRange)));
        yRange   = abs(yAxisMax - yAxisMin);
        yAxisMax = yAxisMax + 0.1*yRange;
        if any(strcmp(MRS_struct.p.target{1}, {'GABA', 'Glx', 'GABAGlx'}))
            yAxisMin = yAxisMin - 0.3*yRange;
        else
            yAxisMin = yAxisMin - 0.1*yRange;
        end
        
    end
    
    axis([0 5 yAxisMin yAxisMax]);
    set(gca,'XDir','reverse','TickDir','out','box','off','XTick',0:5);
    set(get(gca,'YAxis'),'Visible','off');
    xlabel('ppm');
    if MRS_struct.p.HERMES
        title('Difference Spectra and Model Fits');
    else
        title('Difference Spectrum and Model Fit');
    end
    
end



