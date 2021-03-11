function PlotPrePostAlign(MRS_struct, vox, ii, kk)
% Plot pre-/post-alignment spectra
% Updates by MGSaleh 2016, MM 2017-2020

if MRS_struct.p.HERMES
    
    spectraToPlot = zeros(2*length(MRS_struct.p.target), length(MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff));
    count = 0;
    for jj = 1:length(MRS_struct.p.target)
        spectraToPlot((1:2)+count,:) = [MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff(ii,:); ...
                                        MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{jj}).diff_noalign(ii,:)];
        count = count + 2;
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
                if MRS_struct.p.phantom
                    peakRange(jj,:) = MRS_struct.spec.freq <= 4.25 & MRS_struct.spec.freq >= 1.0;
                else
                    peakRange(jj,:) = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.2;
                end
            case 'GSH'
                if ~MRS_struct.p.HERCULES
                    peakRange(jj,:) = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
                else
                    peakRange(jj,:) = MRS_struct.spec.freq <= 4.25 & MRS_struct.spec.freq >= 3;
                end
            case {'Lac','EtOH'}
                peakRange(jj,:) = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
        end
    end
    
    specBaseline  = mean(real(spectraToPlot(:,baseRange)),2);
    spectraToPlot = spectraToPlot - repmat(specBaseline, [1 size(spectraToPlot,2)]);
    
    % Stack spectra
    if MRS_struct.p.phantom
        shift = max(max(real(spectraToPlot(1:2,logical(peakRange(1,:)))))) - ...
                min(min(real(spectraToPlot(1:2,logical(peakRange(1,:))))));
        spectraToPlot(3:4,:) = spectraToPlot(3:4,:) + shift;
        
        yAxisMax = max(real(spectraToPlot(3,logical(peakRange(2,:)))));
        yAxisMin = min(real(spectraToPlot(1,logical(peakRange(1,:)))));
        yRange   = abs(yAxisMax - yAxisMin);
        yAxisMax = yAxisMax + 0.1*yRange;
        yAxisMin = yAxisMin - 0.1*yRange;
    else
        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'}))
            maxGABAGlx           = abs(max(max(real(spectraToPlot(1:2,logical(peakRange(1,:)))),[],2)));
            minGSH               = abs(min(min(real(spectraToPlot(3:4,logical(peakRange(2,:)))),[],2)));
            shift                = max([maxGABAGlx minGSH]) + 0.5*min([maxGABAGlx minGSH]);
            spectraToPlot(3:4,:) = spectraToPlot(3:4,:) + shift;
            
            yAxisMax = max(real(spectraToPlot(3,logical(peakRange(2,:)))));
            yAxisMin = min(real(spectraToPlot(1,logical(peakRange(1,:)))));
            yRange   = abs(yAxisMax - yAxisMin);
            yAxisMax = yAxisMax + 0.1*yRange;
            yAxisMin = yAxisMin - 0.1*yRange;
        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
            shift                = abs(min(real(spectraToPlot(5,logical(peakRange(3,:))))));
            spectraToPlot(5:6,:) = spectraToPlot(5:6,:) + 1.5*shift;
            shift                = abs(max(real(spectraToPlot(1,logical(peakRange(1,:))))));
            spectraToPlot(1:2,:) = spectraToPlot(1:2,:) - 1.5*shift;
            
            yAxisMax = abs(max(real(spectraToPlot(5,logical(peakRange(3,:))))));
            yAxisMax = yAxisMax + 0.2*yAxisMax;
            yAxisMin = max(real(spectraToPlot(1,logical(peakRange(1,:)))));
            yAxisMin = yAxisMin - 4*abs(yAxisMin);
        else
            shift                = abs(min(real(spectraToPlot(5,logical(peakRange(3,:))))));
            spectraToPlot(5:6,:) = spectraToPlot(5:6,:) + 1.5*shift;
        end
    end
    
    count = 0;
    for jj = 1:length(MRS_struct.p.target)
        hold on;
        plot(MRS_struct.spec.freq, real(spectraToPlot(2+count*2,:)), 'Color', 'r');
        plot(MRS_struct.spec.freq, real(spectraToPlot(1+count*2,:)), 'Color', 'b');
        hold off;
        count = count + 1;
    end
    
else
    
    spectraToPlot = [MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff(ii,:); ...
                     MRS_struct.spec.(vox{kk}).(MRS_struct.p.target{1}).diff_noalign(ii,:)];
    
    % Shift baselines to zero
    baseRange = MRS_struct.spec.freq <= 0 & MRS_struct.spec.freq >= -0.5;
    % Some bandwidth-limited acquisitions may not record anything below
    % 0 ppm, in this case get the baseline from the other side of water
    if sum(baseRange) == 0
        baseRange = MRS_struct.spec.freq >= 7 & MRS_struct.spec.freq <= 8;
    end
    switch MRS_struct.p.target{1}
        case {'GABA','Glx','GABAGlx'}
            if MRS_struct.p.phantom
                peakRange = MRS_struct.spec.freq <= 4.25 & MRS_struct.spec.freq >= 1.0;
            else
                peakRange = MRS_struct.spec.freq <= 4.1 & MRS_struct.spec.freq >= 2.26;
            end
        case 'GSH'
            peakRange = MRS_struct.spec.freq <= 3.5 & MRS_struct.spec.freq >= 0.5;
        case {'Lac','EtOH'}
            peakRange = MRS_struct.spec.freq <= 3 & MRS_struct.spec.freq >= 0.5;
    end
    
    specBaseline  = mean(real(spectraToPlot(:,baseRange)),2);
    spectraToPlot = spectraToPlot - repmat(specBaseline, [1 length(spectraToPlot)]);
    
    % Stack spectra
    shift = max(abs(real(spectraToPlot(1,peakRange))));
    spectraToPlot(2,:) = spectraToPlot(2,:) + shift + 0.05*shift;
    
    hold on;
    plot(MRS_struct.spec.freq, real(spectraToPlot(2,:)), 'Color', 'r');
    plot(MRS_struct.spec.freq, real(spectraToPlot(1,:)), 'Color', 'b');
    hold off;
    
    yAxisMax = max(max(real(spectraToPlot(:,peakRange)),[],2));
    yAxisMin = min(min(real(spectraToPlot(:,peakRange)),[],2));
    yRange   = abs(yAxisMax - yAxisMin);
    yAxisMax = yAxisMax + 0.1*yRange;
    if any(strcmp(MRS_struct.p.target{1}, {'GABA','Glx','GABAGlx'}))
        if MRS_struct.p.phantom
            yAxisMin = yAxisMin - 0.15*abs(yAxisMin);
        else
            yAxisMin = yAxisMin - 0.3*yRange;
        end
    else
        yAxisMin = yAxisMin - 0.1*yRange;
    end
    
end

if strcmp(MRS_struct.p.target{1},'EtOH')
    legend({'pre','post'},'EdgeColor',[1 1 1],'Location','northwest');
else
    legend({'pre','post'},'EdgeColor',[1 1 1]);
end
axis([0 5 yAxisMin yAxisMax]);
set(gca,'XDir','reverse','TickDir','out','box','off','XTick',0:5);
set(get(gca,'YAxis'),'Visible','off');
xlabel('ppm');
if MRS_struct.p.HERMES
    title({'Difference Spectra';'(pre- and post-alignment)'});
else
    title({'Difference Spectrum';'(pre- and post-alignment)'});
end



