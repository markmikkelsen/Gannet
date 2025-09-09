function MRS_struct = SpecifyOnOffOrder(MRS_struct)
% Determines the subexperiment order for MEGA-PRESS and HERMES. Borrows
% heavily from osp_onOffClassifyMEGA.m and osp_onOffClassifyHERMES.m from
% Osprey.
%
% Author:
%   Dr. Georg Oeltzschner (Johns Hopkins University, 2018-02-24)
%   goeltzs1@jhmi.edu
%
% History:
%   2018-02-24: First version.
%   2018-11-19: Second version.
%   2020-10-13: Third version.
%   2021-06-28: Fourth version. Adopted code from Osprey.
%   2024-01-22: Fifth version. Reintroduced option to manually select
%               subexperiment order in GannetPreInitialise.m
%
% [1 = ON, 0 = OFF]

if ~isempty(MRS_struct.p.ON_OFF_order)

    if MRS_struct.p.HERMES
        ON_OFF = zeros(2,4);
        for ii = 1:length(MRS_struct.p.ON_OFF_order)
            if strcmpi(MRS_struct.p.ON_OFF_order(ii), 'A')
                ON_OFF(:,ii) = [1 1]';
            elseif strcmpi(MRS_struct.p.ON_OFF_order(ii), 'B')
                ON_OFF(:,ii) = [1 0]';
            elseif strcmpi(MRS_struct.p.ON_OFF_order(ii), 'C')
                ON_OFF(:,ii) = [0 1]';
            elseif strcmpi(MRS_struct.p.ON_OFF_order(ii), 'D')
                ON_OFF(:,ii) = [0 0]';
            end
        end
        MRS_struct.fids.ON_OFF = repmat(ON_OFF, [1 size(MRS_struct.fids.data,2)/4]);
    else
        if strcmpi(MRS_struct.p.ON_OFF_order, 'offfirst')
            MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
        else
            MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
        end
    end

else

    spec      = abs(real(fftshift(fft(MRS_struct.fids.data,[],1),1)));
    freqRange = MRS_struct.p.sw(1) / MRS_struct.p.LarmorFreq(1);
    freq      = (MRS_struct.p.npoints(1) + 1 - (1:MRS_struct.p.npoints(1))) / MRS_struct.p.npoints(1) * freqRange + 4.68 - freqRange/2;
    NAAlim    = freq <= 2.3 & freq >= 1.7;
    waterLim  = freq <= 5 & freq >= 4;
    waterLim2 = freq <= 4.25 & freq >= 3.5;

    switch num2str(length(MRS_struct.p.target))

        case '1'
            if any(strcmp(MRS_struct.p.target, {'GABA', 'Glx', 'GABAGlx'}))
                freqLim = NAAlim;
            elseif strcmp(MRS_struct.p.target, 'GSH')
                freqLim = waterLim;
            elseif any(strcmp(MRS_struct.p.target, {'Lac', 'EtOH'}))
                freqLim = waterLim2;
            end

            specA = mean(spec(freqLim, 1:2:end),2);
            specB = mean(spec(freqLim, 2:2:end),2);

            max_diffAB = max(specA - specB);
            max_diffBA = max(specB - specA);

            if max_diffAB > max_diffBA
                MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
            else
                MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
            end

        case {'2','3'}

            [max_first(1), max_second(1)] = findMaxNAAw(mean(spec(:,1:4:end),2), waterLim, NAAlim);
            [max_first(2), max_second(2)] = findMaxNAAw(mean(spec(:,2:4:end),2), waterLim, NAAlim);
            [max_first(3), max_second(3)] = findMaxNAAw(mean(spec(:,3:4:end),2), waterLim, NAAlim);
            [max_first(4), max_second(4)] = findMaxNAAw(mean(spec(:,4:4:end),2), waterLim, NAAlim);

            % Sort the intensities in ascending order
            [~,order_first]  = sort(max_first);
            [~,order_second] = sort(max_second);

            % Now loop over the subspectra indices (A = 1, B = 2, etc) to determine
            % whether the respective experiments have high or low intensities
            GABA_ON = zeros(1,4);
            GSH_ON  = zeros(1,4);
            EtOH_ON = zeros(1,4);

            for ii = 1:4
                idx_first  = find(order_first == ii);
                idx_second = find(order_second == ii);

                if ismember(idx_second,[3 4])
                    GABA_ON(ii) = 0;
                elseif ismember(idx_second,[1 2])
                    GABA_ON(ii) = 1;
                end

                if ismember(idx_first,[3 4])
                    GSH_ON(ii) = 0;
                elseif ismember(idx_first,[1 2])
                    GSH_ON(ii) = 1;
                end

                if length(MRS_struct.p.target) == 3
                    if (GABA_ON(ii) == 1 && GSH_ON(ii) == 1) || (GABA_ON(ii) == 0 && GSH_ON(ii) == 0)
                        EtOH_ON(ii) = 0;
                    elseif (GABA_ON(ii) == 1 && GSH_ON(ii) == 0) || (GABA_ON(ii) == 0 && GSH_ON(ii) == 1)
                        EtOH_ON(ii) = 1;
                    end
                end
            end

            if length(MRS_struct.p.target) == 2
                MRS_struct.fids.ON_OFF = repmat([GABA_ON; GSH_ON], [1 size(MRS_struct.fids.data,2)/4]);
            elseif length(MRS_struct.p.target) == 3
                MRS_struct.fids.ON_OFF = repmat([EtOH_ON; GABA_ON; GSH_ON], [1 size(MRS_struct.fids.data,2)/4]);
            end

    end

end

end


function [max_w, max_NAA] = findMaxNAAw(spec, waterLim, NAAlim)
% This embedded function finds the maximum intensities of the water
% and NAA signals of an input spectrum.

% Determine maximum absolute signal
max_w   = max([abs(max(spec(waterLim))), abs(min(spec(waterLim)))]);
max_NAA = max([abs(max(spec(NAAlim))), abs(min(spec(NAAlim)))]);

end


% function [max_ins, max_NAA] = findMaxNAACr(in)
% % This embedded function finds the maximum intensities of the creatine and
% % siganls at 3.9 ppm and NAA signals of an input spectrum.
%
% % Determine relevant frequency ranges
% out_ins = op_freqrange(in,3.8,4);
% out_NAA = op_freqrange(in,1.8,2.2);
%
% % Determine maximum absolute signal
% max_ins = max([abs(max(real(out_ins.specs))), abs(min(real(out_ins.specs)))]);
% max_NAA = max([abs(max(real(out_NAA.specs))), abs(min(real(out_NAA.specs)))]);
% end


% switch MRS_struct.p.ON_OFF_order
%
%     case 'onfirst'
%
%         if MRS_struct.p.HERMES
%
%             switch MRS_struct.p.vendor
%
%                 case 'GE'
%                     if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                         if strcmpi(MRS_struct.p.seqorig,'Lythgoe')
%                             % 1=ExpD, 2=ExpB, 3=ExpC, 4=ExpA (MM: 201013)
%                             MRS_struct.fids.ON_OFF = repmat([0 1 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif strcmpi(MRS_struct.p.seqorig,'Noeske')
%                             % 1=ExpB, 2=ExpC, 3=ExpA, 4=ExpD (MM: 201013)
%                             MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         else
%                             % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
%                             MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                         MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]);
%                     elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
%                         % 1=B, 2=C, 3=A, 4=D
%                         MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%                 case {'Philips','Philips_data'}
%                     if ~MRS_struct.p.HERCULES
%                         if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                             % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
%                             MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                             % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
%                             MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     else
%                         % 1=ExpC, 2=ExpD, 3=ExpA, 4=ExpB
%                         MRS_struct.fids.ON_OFF = repmat([0 0 1 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%                 case {'Siemens_twix','Siemens_dicom'}
%                     if ~MRS_struct.p.HERCULES
%                         if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                             % 1=ExpB, 2=ExpD, 3=ExpC, 4=ExpA (MM: 181210)
%                             MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                             % This has not been tested with universal sequence -- 03142018 MGSaleh
%                             MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
%                             % 1=A, 2=D, 3=B, 4=C
%                             MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     else
%                         % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD
%                         MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%             end
%
%         else % MEGA-PRESS
%
%             if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
%                 MRS_struct.fids.ON_OFF = [ones(1,size(MRS_struct.fids.data,2)/2) zeros(1,size(MRS_struct.fids.data,2)/2)];
%             elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
%                 if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
%                     MRS_struct.fids.ON_OFF = repmat([ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) ...
%                         zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))], ...
%                         [1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
%                 else
%                     MRS_struct.fids.ON_OFF = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     %MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]); %This seems to work with the MR1 Philips_data
%                 end
%             else
%                 MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
%             end
%
%         end
%
%     case 'offfirst'
%
%         if MRS_struct.p.HERMES
%
%             switch MRS_struct.p.vendor
%
%                 case 'GE'
%                     if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                         if strcmpi(MRS_struct.p.seqorig,'Lythgoe')
%                             % 1=ExpD, 2=ExpB, 3=ExpC, 4=ExpA (MM: 201013)
%                             MRS_struct.fids.ON_OFF = repmat([0 1 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif strcmpi(MRS_struct.p.seqorig,'Noeske')
%                             % 1=ExpB, 2=ExpC, 3=ExpA, 4=ExpD (MM: 201013)
%                             MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         else
%                             % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
%                             MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                         MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%                 case {'Philips','Philips_data'}
%                     if ~MRS_struct.p.HERCULES
%                         if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                             % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
%                             MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                             MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     else
%                         % 1=ExpC, 2=ExpD, 3=ExpA, 4=ExpB
%                         MRS_struct.fids.ON_OFF = repmat([0 0 1 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%                 case 'Siemens_twix'
%                     if ~MRS_struct.p.HERCULES
%                         if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
%                             % 1=ExpB, 2=ExpD, 3=ExpC, 4=ExpA (MM: 181210)
%                             MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
%                             % This has not been tested with universal sequence -- 03142018 MGSaleh
%                             MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
%                             MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
%                         elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
%                             % 1=?, 2=?, 3=?, 4=?
%                             MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                         end
%                     else
%                         % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD
%                         MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
%                     end
%
%             end
%
%         else % MEGA-PRESS
%
%             if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
%                 MRS_struct.fids.ON_OFF = [zeros(1,size(MRS_struct.fids.data,2)/2) ones(1,size(MRS_struct.fids.data,2)/2)];
%             elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
%                 if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
%                     MRS_struct.fids.ON_OFF = repmat([zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) ...
%                         ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))], ...
%                         [1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
%                 else
%                     %MRS_struct.fids.ON_OFF = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
%                     MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
%                 end
%             else
%                 MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
%             end
%
%         end
%
% end



