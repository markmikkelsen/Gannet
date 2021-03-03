function MRS_struct = SpecifyOnOffOrder(MRS_struct)
%   Determines the subexperiment order for MEGA-PRESS and HERMES.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-02-24)
%       goeltzs1@jhmi.edu
%
%   History:
%       2018-02-24: First version.
%       2018-11-19: Second version.
%       2020-10-13: Third version.

% [1 = ON, 0 = OFF]
switch MRS_struct.p.ON_OFF_order
    
    case 'onfirst'
        
        if MRS_struct.p.HERMES
            
            switch MRS_struct.p.vendor
                
                case 'GE'
                    if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                        if strcmpi(MRS_struct.p.seqorig,'Lythgoe')
                            % 1=ExpD, 2=ExpB, 3=ExpC, 4=ExpA (MM: 201013)
                            MRS_struct.fids.ON_OFF = repmat([0 1 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif strcmpi(MRS_struct.p.seqorig,'Noeske')
                            % 1=ExpB, 2=ExpC, 3=ExpA, 4=ExpD (MM: 201013)
                            MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        else
                            % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
                            MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                        MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]);
                    elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                        % 1=B, 2=C, 3=A, 4=D
                        MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
                case {'Philips','Philips_data'}
                    if ~MRS_struct.p.HERCULES
                        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                            % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
                            MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                            % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
                            MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    else
                        % 1=ExpC, 2=ExpD, 3=ExpA, 4=ExpB
                        MRS_struct.fids.ON_OFF = repmat([0 0 1 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
                case {'Siemens_twix','Siemens_dicom'}
                    if ~MRS_struct.p.HERCULES
                        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                            % 1=ExpB, 2=ExpD, 3=ExpC, 4=ExpA (MM: 181210)
                            MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                            % This has not been tested with universal sequence -- 03142018 MGSaleh
                            MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            % 1=A, 2=D, 3=B, 4=C
                            MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    else
                        % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD
                        MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
            end
            
        else % MEGA-PRESS
            
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [ones(1,size(MRS_struct.fids.data,2)/2) zeros(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
                    MRS_struct.fids.ON_OFF = repmat([ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) ...
                        zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))], ...
                        [1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
                else
                    MRS_struct.fids.ON_OFF = repmat([1 1 0 0], [1 size(MRS_struct.fids.data,2)/4]);
                    %MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]); %This seems to work with the MR1 Philips_data
                end
            else
                MRS_struct.fids.ON_OFF = repmat([1 0], [1 size(MRS_struct.fids.data,2)/2]);
            end
            
        end
        
    case 'offfirst'
        
        if MRS_struct.p.HERMES
            
            switch MRS_struct.p.vendor
                
                case 'GE'
                    if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                        if strcmpi(MRS_struct.p.seqorig,'Lythgoe')
                            % 1=ExpD, 2=ExpB, 3=ExpC, 4=ExpA (MM: 201013)
                            MRS_struct.fids.ON_OFF = repmat([0 1 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif strcmpi(MRS_struct.p.seqorig,'Noeske')
                            % 1=ExpB, 2=ExpC, 3=ExpA, 4=ExpD (MM: 201013)
                            MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        else
                            % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD (MM: 171120)
                            MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                        MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
                case {'Philips','Philips_data'}
                    if ~MRS_struct.p.HERCULES
                        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                            % 1=ExpC, 2=ExpB, 3=ExpA, 4=ExpD
                            MRS_struct.fids.ON_OFF = repmat([0 1 1 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                            MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    else
                        % 1=ExpC, 2=ExpD, 3=ExpA, 4=ExpB
                        MRS_struct.fids.ON_OFF = repmat([0 0 1 1; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
                case 'Siemens_twix'
                    if ~MRS_struct.p.HERCULES
                        if all(ismember(MRS_struct.p.target,{'GABAGlx','GSH'})) || all(ismember(MRS_struct.p.target,{'GABA','GSH'})) || all(ismember(MRS_struct.p.target,{'Glx','GSH'}))
                            % 1=ExpB, 2=ExpD, 3=ExpC, 4=ExpA (MM: 181210)
                            MRS_struct.fids.ON_OFF = repmat([1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        elseif all(ismember(MRS_struct.p.target,{'GSH','Lac'}))
                            % This has not been tested with universal sequence -- 03142018 MGSaleh
                            MRS_struct.fids.ON_OFF  = repmat([0 1 1 0], [1 size(MRS_struct.fids.data,2)/4]); % GSH
                            MRS_struct.fids.ON_OFF2 = repmat([0 1 0 1], [1 size(MRS_struct.fids.data,2)/4]); % Lac
                        elseif all(ismember(MRS_struct.p.target,{'EtOH','GABA','GSH'}))
                            % 1=?, 2=?, 3=?, 4=?
                            MRS_struct.fids.ON_OFF = repmat([1 0 1 0; 1 0 0 1; 0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                        end
                    else
                        % 1=ExpA, 2=ExpB, 3=ExpC, 4=ExpD
                        MRS_struct.fids.ON_OFF = repmat([1 1 0 0; 1 0 1 0], [1 size(MRS_struct.fids.data,2)/4]);
                    end
                    
            end
            
        else % MEGA-PRESS
            
            if strcmpi(MRS_struct.p.vendor,'Philips') && strcmpi(MRS_struct.p.seqorig,'Philips')
                MRS_struct.fids.ON_OFF = [zeros(1,size(MRS_struct.fids.data,2)/2) ones(1,size(MRS_struct.fids.data,2)/2)];
            elseif strcmpi(MRS_struct.p.vendor,'Philips_data') % Hardcode for now. This repmat depends on the value entered for NSA
                if ceil(MRS_struct.p.LarmorFreq(MRS_struct.ii)) > 290
                    MRS_struct.fids.ON_OFF = repmat([zeros(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)) ...
                        ones(1,(MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows))], ...
                        [1 size(MRS_struct.fids.data,2)/((MRS_struct.p.Navg(MRS_struct.ii)/MRS_struct.p.nrows)*2)]); % GABA @ 7T % Changed by MGSaleh  -- 2017
                else
                    %MRS_struct.fids.ON_OFF = repmat([0 0 1 1], [1 size(MRS_struct.fids.data,2)/4]);
                    MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
                end
            else
                MRS_struct.fids.ON_OFF = repmat([0 1], [1 size(MRS_struct.fids.data,2)/2]);
            end
            
        end
        
end

end



