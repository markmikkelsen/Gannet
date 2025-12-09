% script for unfolding dual volume data for MRS on PHILIPS
% vincent - UMCU - apr 2014
%
% the function will generate some folders with output in SDAT/SPAR and txt
%
% input:
% (1) folder example: 'D:\DATA\multiband_mrs\20140320_TB'
%     it expects a 'raw' folder in this folder where the raw data is stored
%     needs for the ref scan: cpx or raw and sin and lab
%     for MRS: raw+sin+lab data
% 
% (2) reconstructed voxel;
% % - 1 for only the original/single voxel location, typically the right one for dual volume scans
% % - 2 for only the second voxel location: MANUALLY PROVIDE OFFSET FOR SECOND VOXEL (mrs_offset2)
% % - 0 or any other number for 2-voxel recon (default)
%
% (3) editing yes/no
% % - 0 for default MRS data, output is summed over averages
% % - 1 for editing data, output is not summed over averages yet
%
%
% Version 1.0: Loads Philips raw PRIAM data and creates a Gannet
% MRS_struct. (GO 11/01/2016)
%
function MRS_struct = PhilipsRawRead(MRS_struct, rawfile, recon_voxel, editing)

ii = MRS_struct.ii;

% Clear previous instances
signalunf = [];

% REMOVE OR PARSE THE THREE INPUT ARGUMENTS
if(nargin == 0)
    disp('ERROR - provide data folder:')
    disp(' recon_dual_volume(''D:\data\spectro\'')')
    return
end
if(nargin == 1)
    %only folder input
    recon_voxel = 1; %recon one voxel
    editing = 0;
end
if(nargin == 2)
    %only folder and voxel input
    editing = 0;
end

if MRS_struct.p.PRIAM == 1 % GO 11/01/2016
    recon_voxel = 3; % reconstruct 2 voxels % GO 11/01/2016
    sense_recon = 1; % do SENSE recon % GO 11/01/2016
else % GO 11/01/2016
    recon_voxel = 1; % reconstruct 1 voxel % GO 11/01/2016
    sense_recon = 0; % do classic coil combination % GO 11/01/2016
end % GO 11/01/2016

%% Setup paths and filenames

% Get path, filename, extension
rawpath = pwd;
raw_file = [rawpath filesep rawfile];
if ~(exist([rawpath filesep 'mfiles'])==7)
    mkdir([rawpath filesep 'mfiles']);
end
% Check whether the input *.raw file exists
if ~(exist(raw_file)==2)
    disp('Unable to locate find the specified *.raw file:')
    disp(raw_file)
    return
end

% Put together the filenames for the reference image and the MRS scan.
clc
disp('=== reference scan ===')
files = dir([rawpath, '/*.cpx']);

if(length(files)>1)
    result = input('select ref scan: ');
    if isempty(result)
        return
    end
else
    result = 1;
end
reffile = files(result).name;


disp('===  MRS scan ===')
% GO 11/01/2016: commented out batch processing, only one scan at a time for now
% files = dir([rawpath,'/*.raw']);
files.name = rawfile;
if(length(files)>1)
    result = input('More than one *.raw file found. Select MRS scans (comma separated): ','s');
else
    result = '1';
end

tmp = regexp(result,'([^ ,:]*)','tokens');
array_MRS_files = str2double(cat(2,tmp{:}));

%%spatial separation RE
vox_sep = -input('What is the separation of the voxels in mm?  ');
vox_ang = input('What is the angulation about FH?  ');


%% some more recon and save options

%take coil info from reference imaging scan, or from water MRS scan
%if from water ref scan, first reconstruct sensitivities of two seperate 
% voxels and save to mfiles/sens1 and sens2
% 0 = reference imaging scan
% 1 = water MRS scan
sens_from_ref_scan = 0; % 1 not implemented yet, GO 11/01/2016

%processing parameters
perform_ecc = 1; %eddy current correction
perform_watersubtraction = 1; %subtraction of scaled water scan for sideband reduction
perform_waterremoval_HSVD = 0; %HSVD water removal

save_files_txt = 1; %store txt file (re,im time domain)
save_files_sdat = 1; %store philips SDAT SPAR file
save_images = 1; %store preprocessing plots (histograms, voxel locations)

%% Loop over list of MRS files selected
for MRS_index = 1:length(array_MRS_files)
    
    % Get name of *.raw file to be processed.
    MRSfile = files(array_MRS_files(MRS_index)).name;

    %% Read reference imaging scan (*.cpx)
    
    % Check whether *.mat files from previous loading are present in the
    % mfiles sub-directory.
    
    % If no such files exists, load the *.cpx file first and save the *.mat files.
    if(~exist([rawpath filesep 'mfiles' filesep 'ref_scan_sense_coil_img.mat'],'file') || ...
            ~exist([rawpath filesep 'mfiles' filesep 'ref_scan_volume_coil_img.mat'],'file'))

        disp('loading reference scan: ' )
        disp([rawpath filesep reffile]);

        img = read_cpx([rawpath filesep reffile]);
        
        % If there is a *.raw file for the reference image, extract the
        % noise from this scan.
        norawfile = 1;
        if exist([rawpath filesep 'ref_scan_sense_coil_img.mat'],'file')
            norawfile = 0;
            [~, info] = read_noise([rawpath filesep reffile(1:end-4) '.raw']);
            noise1 = info.FRC_NOISE_DATA(:,:,1);
            noise2 = info.FRC_NOISE_DATA(:,:,2);
        end
        
        img = permute(img,[5 1 2 4 3]);
        % Image dimensions now are:
        % 1 - number of coils
        % 2-4 - image dimensions x, y, z
        % 5 - number of stacks:
        %   1st stack contains receiver coil images
        %   2nd stack, 1st element is body coil image
        ncoils = size(img,1);
        dimx = size(img,2);
        dimy = size(img,3);
        dimz = size(img,4);

        img2 = double(img(:,:,:,:,2));
        img = double(img(:,:,:,:,1));

        img = img/max(abs(img(:)));
        img2 = img2/max(abs(img2(:)));
        
        % Extract noise levels from outermost slice.
        if norawfile
            disp('no raw file, getting noise from edge of image')
            noise1 = squeeze(img(:,1:30,:,1));
            noise1 = reshape(noise1,[ncoils size(noise1,2)*size(noise1,3)]);
            noise2 = squeeze(img2(:,1:30,:,1));
            noise2 = reshape(noise2,[ncoils size(noise2,2)*size(noise2,3)]);
        end

        save([rawpath filesep 'mfiles' filesep 'ref_scan_sense_coil_img.mat'],'img','ncoils','dimx','dimy','dimz','noise1','noise2');
        img = img2;
        clear img2
        save([rawpath filesep 'mfiles' filesep 'ref_scan_volume_coil_img.mat'],'img','ncoils','dimx','dimy','dimz','noise1','noise2');
        clear img1
        disp('done loading ref scan')
    end

    %% Load the previously established reference image and plot the receiver and volume coil noise
    load([rawpath filesep 'mfiles' filesep 'ref_scan_sense_coil_img.mat'])

    figure(1)
    subplot(2,2,1)
    bar(std(noise1'))
    title(['noise ' num2str(ncoils) 'ch receiver']);
    subplot(2,2,2)
    bar(std(noise2'))
    title(['noise volume coil (ch1)']); %what is in the other 31 channels?

    subplot(2,1,2)
    imagesc(real(noise1*noise1'))
    daspect([1 1 1])

    %% Extract voxel size and orientation data from *.sin files.
    
    % Get dimensions from reference image scan
    refsc_voxel_size = get_sin_voxel_size([rawpath filesep reffile(1:end-4) '.sin']);
    refsc_orientation = get_sin_orientation([rawpath filesep reffile(1:end-4) '.sin']);

    % ap_rl_fh
    refsc_offset = get_sin_voxel_offset([rawpath filesep reffile(1:end-4) '.sin']);
    % Reorder: FH, AP, LR
    refsc_offset1 = [refsc_offset(3); refsc_offset(1); refsc_offset(2)];

    % ap_rl_fh
    mrs_offset = get_sin_voxel_offset([rawpath filesep MRSfile(1:end-4) '.sin']);
    mrs_voxelsize = get_sin_voxel_size([rawpath filesep MRSfile(1:end-4) '.sin']);

    disp('NOT CHECKING ON R/P direction')
    switch refsc_orientation
        case 1 %saggital
            %voxel offsets ref scan:  FH     AP    LR
            img = flip(img,4);

            mrs_offset1 = [ mrs_offset(3); mrs_offset(1); mrs_offset(2)];
            disp('MRS offset 2 for 20140307_vol_KK_hippocampus_7T')
            mrs_offset2 = [ mrs_offset(3) + vox_sep; mrs_offset(1); mrs_offset(2)];

        case 2 %coronal
            img = permute(img,[1 2 4 3]);
            dimx = size(img,2);
            dimy = size(img,3);
            dimz = size(img,4);
            %final order: FH, AP, LR
            refsc_voxel_size = [refsc_voxel_size(1); refsc_voxel_size(3); refsc_voxel_size(2)];
            %final order: FH, AP, LR
            disp('Offset Voxel 1:')
            mrs_offset1 = [ mrs_offset(3); mrs_offset(1); mrs_offset(2)]
            disp('Offset Voxel 2:')
            mrs_offset2 = [ mrs_offset(3); mrs_offset(1)+vox_sep*sin(vox_ang*pi/180); mrs_offset(2)  - vox_sep*cos(vox_ang*pi/180)]

        case 3
            disp('ORIENTATION NOT IMPLEMENTED yet')
            return
    end

    mrs_offset1 = mrs_offset1 - refsc_offset1;
    mrs_offset2 = mrs_offset2 - refsc_offset1;

    lb1 = (mrs_offset1 - mrs_voxelsize/2) ./refsc_voxel_size;
    ub1 = (mrs_offset1 + mrs_voxelsize/2) ./refsc_voxel_size;

    xrange = round([dimx/2 + 0.5 + lb1(1):dimx/2 + 0.5 + ub1(1)]); %FH
    yrange = round([dimy/2 + 0.5 + lb1(2):dimy/2 + 0.5 + ub1(2)]); %AP
    zrange = round([dimz/2 + 0.5 + lb1(3):dimz/2 + 0.5 + ub1(3)]); %LR

    % second voxel: fill in correct values!
    lb2 = (mrs_offset2 - mrs_voxelsize/2) ./refsc_voxel_size;
    ub2 = (mrs_offset2 + mrs_voxelsize/2) ./refsc_voxel_size;
    xrange2 = round([dimx/2 + 0.5 + lb2(1):dimx/2 + 0.5 + ub2(1)]); %FH
    yrange2 = round([dimy/2 + 0.5 + lb2(2):dimy/2 + 0.5 + ub2(2)]); %AP
    zrange2 = round([dimz/2 + 0.5 + lb2(3):dimz/2 + 0.5 + ub2(3)]); %LR

    %% Plot an overlay of the voxel on the body coil image
    SOS = squeeze(sqrt(sum(abs(img).^2,1)));

    gray2 = gray;
    gray2(64,:) = [1 0 0];
    gray2(1,:) = [0 1 0];

    figure(2)
    clf

    SOSdisplay = SOS;
    SOSdisplay(xrange,yrange,zrange) = 2;

    subplot(2,3,1)
    imagesc(squeeze(squeeze(SOSdisplay(round(mean(xrange)),:,:))),[-0.2 1.3])
    title('         TRA          A/L')
    axis off;
    axis square;

    subplot(2,3,2)
    imagesc(flipud(squeeze(SOSdisplay(:,round(mean(yrange)),:))),[-0.2 1.3])
    title('         SAG          H/L')
    axis off;
    axis square;

    subplot(2,3,3)
    imagesc(flipud(squeeze(SOSdisplay(:,:,round(mean(zrange))))),[-0.2 1.3])
    title('         COR          H/P')
    axis off;
    axis square;


    SOSdisplay = SOS;
    SOSdisplay(xrange2,yrange2,zrange2) =-2;
    subplot(2,3,4)
    imagesc(squeeze(SOSdisplay(round(mean(xrange2)),:,:)),[-0.2 1.3])
    title('         TRA          A/L')
    axis off;
    axis square;

    subplot(2,3,5)
    imagesc(flipud(squeeze(SOSdisplay(:,round(mean(yrange2)),:))),[-0.2 1.3])
    title('         SAG          H/L')
    axis off;
    axis square;

    subplot(2,3,6)
    imagesc(flipud(squeeze(SOSdisplay(:,:,round(mean(zrange2))))),[-0.2 1.3])
    title('         COR          H/P')
    axis off;
    axis square;

    colormap(gray2)

    mask = (SOS/max(SOS(:)))>0.04;
    save([rawpath filesep 'mfiles' filesep 'mask.mat'],'mask');
    save([rawpath filesep 'mfiles' filesep 'SOS_sense'],'SOS');

    %% generate sensitivity maps
    sensitivities = zeros(ncoils,dimx,dimy,dimz);
    for c=1:ncoils
        sensitivities(c,:,:,:) = mask.*squeeze(img(c,:,:,:))./SOS;
    end

    %% load MRS data
    if ~exist([rawpath filesep 'mfiles' filesep MRSfile(1:end-4) '.mat'],'file')

        [FID, infow1] = loadLabRaw([rawpath filesep MRSfile]);

        ncoils = double(infow1.dims.nCoils);
        npoints = double(infow1.dims.nKx);
        nmixes = double(infow1.dims.nMixes); % KLC Change nMixes to nRows
        nmeasurements = double(infow1.dims.nMeasurements);
        nrows = infow1.dims.nRows;
        naveragesw1 = nmeasurements*nrows;

        FID = squeeze(FID);
        
        % First dimension: coils (32)
        % Second dimension: data points
        % Third dimension: rows
        % Fourth dimension: Measurements (Water-suppressed ON (1) and OFF (2))
        % Fifth dimension: Mixes (ON and OFF spectra)

        FID = reshape(double(FID),ncoils,npoints,nrows,nmixes,nmeasurements);
        FID = permute(FID,[1 2 3 5 4]);
        FID = reshape(FID,ncoils,npoints,nrows*nmeasurements,nmixes);
        FID = permute(FID,[1 2 4 3]);
        
        figure(14)
        clf
        plot(real(FID(:,1:500,2,1))'); % KLC changed 2 to one in second to last dimension
        ylim([-100 100])
    %     nshift = input('number of zeros: ');
        nshift = 0;
        % nshift = 32;
        disp(['SHIFTING: ' num2str(nshift) ' points'])
        FID = circshift(FID,[0 -nshift 0 0]);

        nwaterfiles = sum(squeeze(sum(sum(abs(FID(:,:,2,:)),1),2))~=0);
        sw = 5e3;

        % downsampling
        tic
        [noversamples, npoints]=get_sin_samples([rawpath filesep MRSfile(1:end-4) '.sin']);

        % FIR filter
        clear C
        C(1)= -2.74622659936243e-2;
        C(2)= -3.16256052372452e-2;
        C(3)=  6.00081666037503e-3;
        C(4)=  9.10571641138698e-2;
        C(5)= 1.94603313644247e-1;
        C(6)=  2.66010793957549e-1;
        C(7)=  2.66010793957549e-1;
        C(8)=  1.94603313644247e-1;
        C(9)=  9.10571641138698e-2;
        C(10)=  6.00081666037503e-3;
        C(11)= -3.16256052372452e-2;
        C(12)= -2.74622659936243e-2;

        N_tabs = length(C);
        N_delay = floor(N_tabs/2);
        R_fir = 4;
        R_hdf = 1. * noversamples / npoints / 4;
        % naveragesw1=nrows*nmeasurements;
        % two step undersampling, first R_hdf undersampling
        H_ref_out = squeeze(sum(reshape(FID,[ncoils,R_hdf,noversamples/R_hdf nmixes naveragesw1]),2));

        % second, FIR filter with coefficients from C
        S_ref_out = zeros(ncoils,npoints,nmixes,naveragesw1);

        n = 0:npoints-1;
        n1 = n*R_fir+1+N_delay;
        nmax = noversamples/R_hdf-1;
        for k=1:N_tabs
            n2 = n1-k;
            idx = find(n2>=0 & n2<nmax);
            S_ref_out(:,n(idx)+1,:,:) = S_ref_out(:,n(idx)+1,:,:) + C(k)*H_ref_out(:,n2(idx)+1,:,:);
        end
        FID = S_ref_out;
        disp(['downsampling: ' num2str(toc) ' sec'])

        save([rawpath filesep 'mfiles' filesep MRSfile(1:end-4) '.mat'],'FID','sw','naveragesw1','nwaterfiles','npoints','ncoils');
    else
        load([rawpath filesep 'mfiles' filesep MRSfile(1:end-4) '.mat'])
    end

    %% Calculate MRS sensitivities from first 20 points of water signal (alternative to SENSE reference image)
    % sens1 = mean(FID(:,1:20,2),2);
    % save([folder filesep 'mfiles' filesep 'sens1.mat'],'sens1');
    % %% save MRS sensitivities
    % sens2 = mean(FID(:,1:20,2),2);
    % save([folder filesep 'mfiles' filesep 'sens2.mat'],'sens2');



    %% generate SENSE matrix

    sens1 = squeeze(mean(mean(mean(sensitivities(:,xrange,yrange,zrange),2),3),4));
    sens2 = squeeze(mean(mean(mean(sensitivities(:,xrange2,yrange2,zrange2),2),3),4));
    %get sensitivities
    if ~sens_from_ref_scan

        if exist([rawpath filesep 'mfiles' filesep 'sens1.mat'],'file')
            load([rawpath filesep 'mfiles' filesep 'sens1.mat'],'sens1');
        else
            disp(' !!!    NO SENS1 MATRIX FOUND   !!!')
            disp(' !!!      using ref-scan        !!!')
        end
        if exist([rawpath filesep 'mfiles' filesep 'sens2.mat'],'file')
            load([rawpath filesep 'mfiles' filesep 'sens2.mat'],'sens2');
        else
            disp(' !!!    NO SENS2 MATRIX FOUND   !!!')
            disp(' !!!      using ref-scan        !!!')
        end

        ix = round(mean(xrange));
        iy = round(mean(yrange));
        iz = round(mean(zrange));
        ix2 = round(mean(xrange2));
        iy2 = round(mean(yrange2));
        iz2 = round(mean(zrange2));

        signal = sens1;
        signal = signal./sqrt(sum(abs(signal).^2));
        tmpimage = squeeze(sum(bsxfun(@times, sensitivities, conj(signal)),1));

        disp( ['orig. voxel positions; x1=' num2str(ix) '; y1=' num2str(iy) '; z1=' num2str(iz) ';'])
        [max_val, position] = max(abs(tmpimage(:)));
        [tx,ty,tz] = ind2sub(size(tmpimage),position);
        disp( ['opt. voxel positions; x1=' num2str(tx) '; y1=' num2str(ty) '; z1=' num2str(tz) ';'])
        disp(['similarity: ' num2str(max_val) ' vs ' num2str(abs(tmpimage(ix,iy,iz)))])

        figure(17)
        clf
        jet2 = jet;
        jet2(1,:) = [0;0;0];
        colormap(jet2)

        minpl = 0;
        maxpl = 1;

        subplot(2,3,1)
        imagesc(abs(squeeze(tmpimage(ix,:,:))),[minpl maxpl])
        line([0 100],[iy iy],'Color',[1 1 1])
        line([iz iz],[0 100],'Color',[1 1 1])
        ylabel('y');ylabel('FH')
        xlabel('z');xlabel('LR')
        axis square;

        title('MPR')
        subplot(2,3,2)
        imagesc(abs(squeeze(tmpimage(:,iy,:))),[minpl maxpl])
        line([0 100],[ix ix],'Color',[1 1 1])
        line([iz iz],[0 100],'Color',[1 1 1])
        ylabel('x');ylabel('AP')
        xlabel('z');xlabel('LR')
        axis square;
        set(gca,'YDir','normal');
        title('MPR')
        subplot(2,3,3)
        imagesc(abs(squeeze(tmpimage(:,:,iz))),[minpl maxpl])
        line([0 100],[ix ix],'Color',[1 1 1])
        line([iy iy],[0 100],'Color',[1 1 1])
        ylabel('x');ylabel('AP')
        xlabel('y');xlabel('FH')
        axis square;
        set(gca,'YDir','normal');
        title('MPR')

        signal = sens2;
        signal = signal./sqrt(sum(abs(signal).^2));
        tmpimage = squeeze(sum(bsxfun(@times, sensitivities, conj(signal)),1));

        disp( ['orig. voxel positions; x1=' num2str(ix2) '; y1=' num2str(iy2) '; z1=' num2str(iz2) ';'])
        [max_val, position] = max(abs(tmpimage(:)));
        [tx,ty,tz] = ind2sub(size(tmpimage),position);
        disp( ['opt. voxel positions; x2=' num2str(tx) '; y2=' num2str(ty) '; z2=' num2str(tz) ';'])
        disp(['similarity: ' num2str(max_val) ' vs ' num2str(abs(tmpimage(ix2,iy2,iz2)))])
        disp(' ')

        subplot(2,3,4)
        imagesc(abs(squeeze(tmpimage(ix2,:,:))),[minpl maxpl])
        line([0 100],[iy2 iy2],'Color',[1 1 1])
        line([iz2 iz2],[0 100],'Color',[1 1 1])
        axis square;
        ylabel('y')
        xlabel('z')

        title('MPR')
        subplot(2,3,5)
        imagesc(abs(squeeze(tmpimage(:,iy2,:))),[minpl maxpl])
        line([0 100],[ix2 ix2],'Color',[1 1 1])
        line([iz2 iz2],[0 100],'Color',[1 1 1])
        axis square;
        ylabel('x')
        xlabel('z')
        set(gca,'YDir','normal');

        title('MPR')
        subplot(2,3,6)
        imagesc(abs(squeeze(tmpimage(:,:,iz2))),[minpl maxpl])
        line([0 100],[ix2 ix2],'Color',[1 1 1])
        line([iy2 iy2],[0 100],'Color',[1 1 1])
        axis square;
        ylabel('x')
        xlabel('y')
        set(gca,'YDir','normal');
        title('MPR')

    end

    %normalize
    sens1 = sens1./squeeze(sqrt(sum(abs(sens1).^2,1)));
    sens2 = sens2./squeeze(sqrt(sum(abs(sens2).^2,1)));

    % form S matrix
    if recon_voxel == 1
        S = [sens1];
    elseif recon_voxel == 2
        S = [ sens2];
    else
        S = [ sens1 ...
              sens2 ...
        ];
    end

    PSY = noise1*noise1';
    U = inv(S'*inv(PSY)*S)*S'*inv(PSY);

    ga = inv(S'*inv(PSY)*S);
    gb = S'*inv(PSY)*S;
    g = sqrt(ga.*gb);

    disp(['gfactors: ' num2str(diag(real(g))')])

    %% perform SENSE unfolding 
    % m=1 %1 for metab, 2 for water ref

    disp('sense unfolding...');

    % clear signalunf signalres
    for m=1:2
        for a =1:naveragesw1
            signal = FID(:,:,m,a);
            signalunf(:,:,m,a) = U*signal;
            signalres(:,:,m,a) = signal - S*signalunf(:,:,m);
        end
    end

    nspec = size(signalunf,1);

    %% water subtraction
    searchwindow = [-200 200]; %Hz
    searchwindow = (searchwindow/sw+0.5)*npoints;

    waterfid = mean(signalunf(:,:,2,1:nwaterfiles),4);

    if perform_watersubtraction
        % waterremoval - subtract scaled water scan
        for t=1:nspec
            %
            nruns = 0;

            apo = 5;
            t1 = [0:npoints-1]/sw;
            filt = exp(-t1*apo);


            while(nruns < 5)
                nruns = nruns+1;


                specw = fftshift(fft(waterfid(t,:,:,:).*filt,[],2),2);
                specm = fftshift(fft(mean(signalunf(t,:,1,:),4).*filt,[],2),2);
                [a b] = max(abs(specw));
                [a2 b2] = max(abs(specm));
                corr = specm(b2)./specw(b);

    %             figure(7)
    %             clf

    %             plot([1:2048],circshift(real(specw*corr),[0 b2-b]),[1:2048],real(mean(specm,4)))

                if(b2>max(searchwindow) || b2<min(searchwindow))
    %                 disp([num2str(nruns) ' out of range, stopping'])
    %                 pause(1)
                    continue
                else
    %                 disp([num2str(nruns) ' subtracting'])
    %                 pause(2)
                end

                %
                for a = 1:naveragesw1
                    signalunf(t,:,1,a) = signalunf(t,:,1,a) - corr.*ifft(ifftshift(circshift(specw,[0 b2-b]),2),[],2).*filt;
                end
            end
        end
    end

    %% eddy current correction
    if perform_ecc
    %     w_unf = mean(signalunf(:,:,2,1:nwaterfiles),4);
        w_unf = waterfid;
        wcorr = abs(w_unf)./w_unf;
        wcorr(isnan(wcorr)) = 1;
        for a=1:naveragesw1
            signalunf(:,:,1,a) = signalunf(:,:,1,a).*wcorr;
        end
        for a=1:nwaterfiles
            signalunf(:,:,2,a) = signalunf(:,:,2,a).*wcorr;
        end
    end

    %% water removal
    if perform_waterremoval_HSVD
        % waterremoval
        for t=1:nspec
            for a=1:naveragesw1
              signalunf(t,:,1,a) = waterremovalSVD(signalunf(t,:,1,a).', sw/1000, 32, -0.12, 0.12, 0).';
            end
        end
    end

    %% save spectra to disk 
    save([rawpath filesep 'mfiles' filesep MRSfile(1:end-4) '.mat'],'signalunf','-append'); % GO 02/18/16 save signal output to mat file (easier for GannetLoad to process)

    if ~(exist([rawpath filesep 'rec_spectra_txt'])==7)
        mkdir([rawpath filesep 'rec_spectra_txt']);
    end
    if ~(exist([rawpath filesep 'rec_SDATSPAR'])==7)
        mkdir([rawpath filesep 'rec_SDATSPAR']);
    end

    for m=1:2

        if m==2 ;%water
            saveaverages = nwaterfiles;
        else %metabolites
            saveaverages = naveragesw1;
        end

        for s=1:nspec
    %         for a =1:naveragesw1

                %water or metabolite file
                if(m==2) 
                    watermetstr = 'ref';
                else
                    watermetstr = 'act';
                end

                %voxel number 1 for original, 2 for end one or no number for only one
                if (recon_voxel == 2) || (recon_voxel ==1)
                    savefile = [MRSfile(1:end-4) '_unf_' watermetstr]; %reconstructing only 1 voxel
                else
                    savefile = [MRSfile(1:end-4) '_unf_' watermetstr '_' num2str(s)];
                end

                signalsave = signalunf(s,:,m,1:saveaverages);

                specsave = fftshift(fft(signalsave,[],2),2);

                if(m==1)
                    if editing
                        NAA = max(abs(mean(specsave(1,550:680,:,:),4)));
                    else
                        NAA = max(abs(mean(specsave(1,550:680,:,1:2:end),4)));
                    end
                    disp(['NAA peak intensity: ' num2str(NAA)])
                    tmp = reshape(mean(signalsave,4),256,npoints/256);
                    tmp = std(tmp,[],1);
                    noise = min(real(tmp(:)));
                    disp(['noise (time domain): ' num2str(noise)])

                    noiserec = U*noise1;
                    disp(['noise (FRCnoise): ' num2str(std(noiserec(s,:)))])

                else
                    figure(5)
                    water = max(abs(sum(specsave,4)));
                    f = polyfit([2:15],abs(sum(signalsave(1,2:15,1,:),4)),1);
                    plot([1:15],abs(sum(signalsave(1,1:15,1,:),4)),[1:15],polyval(f,[1:15]));
                    water = polyval(f,1);

                    disp(['water peak intensity (t): ' num2str(water)])
                end


                if ~editing %just store the mean
                    saveaverages = 1;
                    signalsave = mean(signalsave,4);
                    specsave = mean(specsave,4);
                end

                %save txt file
                if save_files_txt
                    for n=1:saveaverages

    %                     if saveaverages ~= 1
                            savefile1 = [savefile '_' num2str(n,'%03d')];
    %                     else
    %                         savefile1 = savefile;
    %                     end

                        disp(['saving: ' savefile1 '.txt'])

                        fid = fopen([rawpath filesep 'rec_spectra_txt' filesep savefile1 '.txt'],'w');
                        if fid == -1
                            disp(['Cannot open file: ' savefile1])
                            continue
                        end
                        for p=1:npoints
                            fprintf(fid,'%f \t %f\r\n',real(signalsave(1,p,1,n)),imag(signalsave(1,p,1,n)));
                        end
                        fclose(fid);
                    end
                end


                if save_files_sdat
                    %save sdatspar
                    fname_p = 'example.SPAR'; %check if this is in the right location
                    fname_out = [rawpath filesep 'rec_SDATSPAR' filesep savefile '.sdat'];
                    disp(['saving: ' fname_out])
                    fname_out_spar = [rawpath filesep 'rec_SDATSPAR' filesep savefile '.spar'];
                    write_sdat_spar(signalsave, fname_p, fname_out, fname_out_spar, saveaverages, npoints)
                end
    %         end
        end
    end

    %% plot unfolded voxels
    % if editing
    %     ncols = 3;
    % else
        ncols = 2;
    % end

    disp(['naverages = ' num2str(naveragesw1)])
    disp(['nwaterfiles = ' num2str(nwaterfiles)])
    % averages
    m_unfdisp = signalunf(:,:,1,:);
    w_unfdisp = signalunf(:,:,2,1:nwaterfiles);

    figure(4);
    subplot(2,1,1)
    imagesc(abs(squeeze(m_unfdisp(1,:,:,:))))
    subplot(2,1,2)
    imagesc(abs(squeeze(w_unfdisp(1,:,:,:))))


    m_unfdisp = mean(m_unfdisp(:,:,:,:),4);
    w_unfdisp = mean(w_unfdisp(:,:,:,:),4);


    % noise = std(signalunfdisp(:,1000:1500,1),[],2).';
    % disp(['noise = ' num2str(noise)])
    % for t=1:nspec
    %     signalunfdisp(t,:,:) = signalunfdisp(t,:,:)./noise(t);
    % end   
    % for t=1:nspec
    %     signalunfdisp(t,:,1) = signalunfdisp(t,:,1)/max(abs(signalunfdisp(t,:,1)));
    % end
    zff = 2;
    apo = 1;
    t = [0:npoints-1]/sw;
    fs = 7;
    ppm = ([0:npoints*zff-1]/(npoints*zff-1)-0.5)*sw / (298/7*fs)+4.68;

    filt = exp(-t*apo);
    filt = repmat(filt,[nspec 1]);
    % specunfodd = fftshift(fft(signalunfdispodd.*filt,[],2),2);
    % specunfeven = fftshift(fft(signalunfdispeven.*filt,[],2),2);
    m_specunf = fftshift(fft(m_unfdisp.*filt,[npoints*zff],2),2);
    w_specunf = fftshift(fft(w_unfdisp.*filt,[npoints*zff],2),2);

    % angle(signalunfdisp(1,1))/pi*180

    ph(1) = 0;
    ph(2) = 0;

    figure(3)
    clf

    xmin = 4;
    % xmax = 5.36;
    xmax = 5.5;

    % plot water
    % axis1 = [xmin xmax 1.1*min([real(w_specunf(:)); imag(w_specunf(:))]) 1.1*max([real(w_specunf(:)); imag(w_specunf(:))])];
    axis1 = [xmin xmax 1.1*min(real(w_specunf(:))) 1.1*max(real(w_specunf(:)))];
    % axis1 = [xmin xmax -1.5e5 4e5];

    for t=1:nspec
        subplot(ncols,nspec,t)

    %     plot(ppm,imag(exp(-1i*ph(t))*w_specunf(t,:)),'g')
    %     hold on
        plot(ppm,real(exp(-1i*ph(t))*w_specunf(t,:)))
    %     hold off
        set(gca,'XDir','reverse');
        axis(axis1);
        title(['water spectrum, voxel: ' num2str(t)])

    %     disp(['water signal intensity: ' num2str(max(abs(w_specunf(t,:))))])
    end

    xmin = 0;
    xmax = 4.2;

    % plot metabolites
    axis1 = [xmin xmax 1.1*min(real(m_specunf(:))) 1.1*max(real(m_specunf(:)))];

    for t=1:nspec
        subplot(ncols,nspec,t+nspec)
        plot(ppm,real(exp(-1i*ph(t))*m_specunf(t,:)))
        set(gca,'XDir','reverse');
        axis(axis1)
        title(['metabolite spectrum, voxel: ' num2str(t)])
    end

    %% store images 
    if save_images
        if ~exist([rawpath filesep 'images'],'dir')
            mkdir([rawpath filesep 'images']);
        end

        % store figures
        filename = [rawpath filesep 'images' filesep MRSfile(1:end-4)];
        figure(1)
        file = [filename '_noise'];
        print('-vector','-r300','-dpng',file);
        print('-vector','-depsc2',file);


        figure(2)
        file = [filename '_location'];
        print('-vector','-r300','-dpng',file);
        print('-dpsc2', '-noui', '-vector', file);
        print('-vector','-depsc2',file);

        figure(3)
        if(recon_voxel == 2) || (recon_voxel ==1)
            file = filename;
        else
            file = [filename '_unfolded'];
        end
        print('-vector','-r300','-dpng',file);
        print('-vector','-depsc2',file);
    end
    
end
disp(' ')
disp('recon_dual_volume: finished')

% FIDs may be in the wrong order depending on NMixes, correct that here for now, fix later %
% GO 11/03/2016
signalunf(:,:,:,2:2:end) = -signalunf(:,:,:,2:2:end);
a = signalunf(:,:,:,1:end/2);
b = signalunf(:,:,:,1+end/2:end);
c = zeros(size(signalunf));
c(:,:,:,1:2:end) = a;
c(:,:,:,2:2:end) = b;
signalunf = c;
clear a b c;
a = zeros(size(signalunf));
a(:,:,:,1:4:end) = signalunf(:,:,:,1:4:end);
a(:,:,:,2:4:end) = signalunf(:,:,:,3:4:end);
a(:,:,:,3:4:end) = signalunf(:,:,:,2:4:end);
a(:,:,:,4:4:end) = signalunf(:,:,:,4:4:end);
signalunf = a;
clear a;

% Save all relevant data/information to MRS_struct % GO 11/01/2016
MRS_struct.p.NVoxels                = size(signalunf,1);
MRS_struct.p.npoints(ii)            = npoints; % GO 11/01/2016
MRS_struct.p.nrows(ii)              = size(signalunf,4); % GO 11/01/2016
MRS_struct.p.ncoils                 = ncoils; % GO 11/01/2016
MRS_struct.p.Navg(ii)               = size(signalunf,4); % GO 11/01/2016
MRS_struct.p.Nwateravg(ii)          = nwaterfiles; % GO 11/01/2016
MRS_struct.p.voxdim(ii,:)           = [mrs_voxelsize(1) mrs_voxelsize(2) mrs_voxelsize(3)]; %AP, RL, FH - preliminary, TEST! % GO 11/01/2016
MRS_struct.p.voxoff(ii,:)           = [mrs_offset(1) mrs_offset(2) mrs_offset(3)]; %AP, RL, FH - preliminary, TEST! % GO 11/01/2016
MRS_struct.p.voxang(ii,:)           = vox_ang; % voxel angulation (1 dimension only so far) % GO 11/01/2016
MRS_struct.p.TR(ii)                 = get_sin_TR([rawpath filesep MRSfile(1:end-4) '.sin']);% GO 11/01/2016
MRS_struct.p.TE(ii)                 = get_sin_TE([rawpath filesep MRSfile(1:end-4) '.sin']);% GO 11/01/2016
MRS_struct.p.LarmorFreq(ii)         = 127; % Need to get that from somewhere! GO 11/01/2016
MRS_struct.p.sw(ii)                 = 2e3; % Need to parse that from somewhere! GO 11/01/2016
MRS_struct.multivoxel.sensitivities = sensitivities; % GO 11/01/2016
MRS_struct.multivoxel.voxsep        = vox_sep; % voxel separation (1 dimension only so far) % GO 11/01/2016

% save all unfolded signals to MRS_struct
MRS_struct.multivoxel.allsignals = signalunf; % GO 11/01/2016

% work up water data to plug back into GannetLoad % GO 11/02/2016
MRS_struct.fids.data_water = conj(squeeze(MRS_struct.multivoxel.allsignals(:,:,2,1:MRS_struct.p.Nwateravg))); % select non-zero water spectra % GO 11/02/2016
MRS_struct.out.phase_water = conj(MRS_struct.fids.data_water(1))./abs(MRS_struct.fids.data_water(1));

save([pwd filesep 'MRS_struct.mat'],'MRS_struct')

end
%% END OF RECONSTRUCTION CODE





%% ADDITIONAL FUNCTIONS FOR DATA PROCESSING

%%display file types
function a = display_file_names(folder,filetype)

if ~(exist([folder filesep 'mfiles'])==7)
    mkdir([folder filesep 'mfiles']);
end

a = dir([folder filesep 'raw\*' filetype]);
for n=1:length(a)
    disp(['(' num2str(n) ') file = ''' a(n).name(1:end-4) ''';'])
end

end

%%read sin file voxel size
function [voxel_sizes]=get_sin_voxel_size(filename)

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    voxel_sizes = -1;
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'voxel_sizes') ~= 0)
        if line(12:22)=='voxel_sizes'
            voxel_sizes = sscanf(line, '%*s %*s %*s %*s %*s %f %f %f');
        end
    end
end
fclose(sinfid);
end

%%read sin file voxel offset
function [voxel_offset]=get_sin_voxel_offset(filename)

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    voxel_sizes = -1;
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'loc_ap_rl_fh_offcentres') ~= 0)
        if line(12:34)=='loc_ap_rl_fh_offcentres'
            voxel_offset = sscanf(line, '%*s %*s %*s %*s %*s %f %f %f');
        end
    end
%     if(strfind(line,'loc_ap_rl_fh_offcentr_incrs') ~= 0)
%         if line(12:38)=='loc_ap_rl_fh_offcentr_incrs'
%             voxel_offset = sscanf(line, '%*s %*s %*s %*s %*s %f %f %f');
%         end
%     end
    
end
fclose(sinfid);
end

%%read sin file voxel offset
function [orientation]=get_sin_orientation(filename)

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    voxel_sizes = -1;
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'slice_orientations') ~= 0)
        if line(12:29)=='slice_orientations'
            orientation = sscanf(line, '%*s %*s %*s %*s %*s %d');
        end
    end
end
fclose(sinfid);
end

% read sin TR % GO 11/01/2016
function TR=get_sin_TR(filename)

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'repetition_times') ~= 0)
        if line(12:27)=='repetition_times'
            TR = sscanf(line, '%*s %*s %*s %*s %*s %d');
        end
    end
end
fclose(sinfid);
end

% read sin TE % GO 11/01/2016
function TE=get_sin_TE(filename)

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'echo_times') ~= 0)
        if line(12:21)=='echo_times'
            TE = sscanf(line, '%*s %*s %*s %*s %*s %d');
        end
    end
end
fclose(sinfid);
end

%%read sin file samples
function [cMRS_samples, cMRS_samples_ds]=get_sin_samples(filename)
cMRS_samples = -1;
cMRS_samples_ds = -1;

toks = regexpi(filename,'^(.*?)(\.sin|\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
sinfilename = sprintf('%s.sin',prefix);
sinfid=fopen(sinfilename);
if(sinfid == -1)
    disp(['file not found: ' sinfilename])
    return
end
while ~feof(sinfid)
    line=fgetl(sinfid);
    if(strfind(line,'max_dr_samples') ~= 0)
        if line(12:25)=='max_dr_samples'
            cMRS_samples = sscanf(line, '%*s %*s %*s %*s %*s %d');
        end
    end
    if(strfind(line,'dc_freq_domain_length') ~= 0)
        if line(12:32)=='dc_freq_domain_length'
            cMRS_samples_ds = sscanf(line, '%*s %*s %*s %*s %*s %d');
        end
    end    
end
fclose(sinfid);
end

%%read cpx
function data = read_cpx(file, border, flip_img, kspace, read_params, compression_parameter)

%--------------------------------------------------------------------------
% data = read_cpx(file, border, flip_img)
%
% read_cpx: Function for reading the whole cpx-file and writes the data in
%           a predefined 9 dimensional array
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               border      0 returns the original sized image
%                           1 makes a square image
%               flip_img    0 no flip of the image
%                           1 flips the image (like the scanner does)
%               read_params optional input. Specifies the images to be
%                           read. read_params is a struct created
%                           by the function create_read_param_struct
%
% Output:       data        A 9 dimensional array which consists all the
%                           data in the cpx file. The array has the
%                           following structure:
%                           Column: 1:resolution y, 2:resolution x 3:Stack,
%                                   4:Slice, 5:Coil, 6:Heart phase, 7:Echo,
%                                   8:Dynamics, 9:Segments, 10:Segments2
%
% Notice! For scans with body coil and surface coil in different stacks,
% e.g. SENSE reference scan, the body coil is numbered as coil no. 1 in
% stack no. 2, the surface coils are numbered as coils 2,3,4... in stack no. 1!
%------------------------------------------------------------------------


switch nargin
    case 5
        compression_parameter = [];
    case 4
        read_params = create_read_param_struct(file);
        compression_parameter = [];
    case 3
        read_params = create_read_param_struct(file);
        kspace = 0;
        compression_parameter = [];
    case 2
        flip_img = 0;
        kspace = 0;
        read_params = create_read_param_struct(file);
        compression_parameter = [];
    case 1
        flip_img = 0;
        border = 0;
        kspace = 0;
        read_params = create_read_param_struct(file);
        compression_parameter = [];
end

%Reads the header of the file
header = read_cpx_header(file,'no');
[rows,columns] = size(header);

% Calculates the number of slices, coils, etc...
stacks   = length(read_params.loca);
slices   = length(read_params.slice);
coils    = length(read_params.coil);
hps      = length(read_params.phase);
echos    = length(read_params.echo);
dynamics = length(read_params.dyn);
segments = length(read_params.seg);
segments2= length(read_params.seg2);


res_x       = header(1,9);
res_y       = header(1,10);
compression = header(1,11);
flip        = header(1,12);

offset_table_cpx=create_offset_table(header);

% defines the array for the output
if border
    res = max([res_x, res_y]);
    res_x = res;
    res_y = res;
end

if ~isempty(compression_parameter)
    data = zeros(res_x, res_y, stacks, slices, compression_parameter{1}, hps, echos, dynamics, segments, segments2,'single');
    data3 = zeros(res_x, res_y,coils);
else
    data = zeros(res_x, res_y, stacks, slices, coils, hps, echos, dynamics, segments, segments2,'single');
    data3 = zeros(res_x, res_y,coils);
end

% Define a waitbar
h = waitbar(0, 'Loading file...');
set(h,'Units','pixels')
scnsize = get(0,'ScreenSize');
set(h,'Position',[floor(scnsize(3)/2)-160,floor(scnsize(4)/2)-30,360,75]);

fid = fopen(file);

% Runs through all images in the file, reads them and writes it in the
% correct position in the output array "data"
i = 1;
total_loops = 1;
for loop = 1:2
    for st = 1:stacks
        for sl = 1:slices
            for se2 = 1:segments2
                for ph = 1:hps
                    for ec = 1:echos
                        for dy = 1:dynamics
                            for se = 1:segments
                                for co = 1:coils
                                    offset = offset_table_cpx(read_params.loca(st),read_params.slice(sl),read_params.coil(co),read_params.phase(ph),read_params.echo(ec),read_params.dyn(dy),read_params.seg(se),read_params.seg2(se2));
                                    if offset >=0
                                        if loop == 2
                                            image = read_cpx_image(file, offset, border, flip_img);
                                            if kspace
                                                image = fftshift(fft2(fftshift(image)));
                                            end
                                            data3(:,:,co) = image;
                                            waitbar(i/total_loops,h)
                                            i = i+1;
                                        else
                                            total_loops = total_loops +1;
                                        end
                                    end
                                end
                                if ~isempty(compression_parameter)
                                    %                                     data_temp= squeeze(combine_data(reshape(data3,size(data3,1),size(data3,2), 1,size(data3,3)),compression_parameter{2}));
                                    data(:,:,st,sl,:,ph,ec,dy,se, se2) = reshape(combine_data_gui(reshape(data3,size(data3,1),size(data3,2), 1,size(data3,3)),compression_parameter{2}),size(data3,1),size(data3,2),1,1,compression_parameter{1});
                                else
                                    data(:,:,st,sl,:,ph,ec,dy,se, se2) = reshape(data3,size(data3,1),size(data3,2),1,1,size(data3,3));
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

close(h);
fclose all;

end

%%create read param struct
function [v,raw_params] = create_read_param_struct(file)

dotind = strfind(file,'.');
ending = lower(file(dotind(end)+1:end));

switch ending
    case 'rec'
        parfile = [file(1:dotind),'par'];
        par = parread(parfile);
        v.slice = unique(par.ImageInformation.SliceNumber);
        v.echo = unique(par.ImageInformation.EchoNumber);
        v.dyn = unique(par.ImageInformation.DynamicScanNumber);
        v.phase =  unique(par.ImageInformation.CardiacPhaseNumber);
        v.typ = unique(par.ImageInformation.ImageTypeMr);
        raw_params = par;
    case 'par'
        par = parread(file);
        v.slice = unique(par.ImageInformation.SliceNumber);
        v.echo = unique(par.ImageInformation.EchoNumber);
        v.dyn = unique(par.ImageInformation.DynamicScanNumber);
        v.phase =  unique(par.ImageInformation.CardiacPhaseNumber);
        v.typ = unique(par.ImageInformation.ImageTypeMr);
        raw_params = par;
    case 'cpx'
        header = read_cpx_header(file,'no');
        v.loca = unique(header(:,1))+1;
        v.slice = unique(header(:,2))+1;
        v.coil = unique(header(:,3))+1;
        v.phase = unique(header(:,4))+1;
        v.echo = unique(header(:,5))+1;
        v.dyn = unique(header(:,6))+1;
        v.seg = unique(header(:,7))+1;
        v.seg2 = unique(header(:,18))+1;
        raw_params = header;
    case 'data'
        t = 'TEHROA';
        listfile = [file(1:dotind),'list'];
        list = listread(listfile);
        typ = unique(list.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = strfind(typ(i),t);
        end
        v.typ = sort(numtyp);
        v.mix = unique(list.Index.mix)+1;
        v.dyn = unique(list.Index.dyn)+1;
        v.phase = unique(list.Index.card)+1;
        v.echo = unique(list.Index.echo)+1;
        v.loca = unique(list.Index.loca)+1;
        v.coil = unique(list.Index.chan)+1;
        v.seg = unique(list.Index.extr1)+1;
        v.seg2 = unique(list.Index.extr2)+1;
        v.ky = unique(list.Index.ky)+1;
        v.slice = unique(list.Index.kz)+1;
        v.aver = unique(list.Index.aver)+1;
        raw_params = list;
    case 'list'
        t = 'TEHROA';
        list = listread(file);
        typ = unique(list.Index.typ(:,2));
        for i = 1:length(typ)
            numtyp(i) = strfind(typ(i),t); %#ok<*AGROW> 
        end
        v.typ = sort(numtyp);
        v.mix = unique(list.Index.mix)+1;
        v.dyn = unique(list.Index.dyn)+1;
        v.phase = unique(list.Index.card)+1;
        v.echo = unique(list.Index.echo)+1;
        v.loca = unique(list.Index.loca)+1;
        v.coil = unique(list.Index.chan)+1;
        v.seg = unique(list.Index.extr1)+1;
        v.seg2 = unique(list.Index.extr2)+1;
        v.ky = unique(list.Index.ky)+1;
        v.slice = unique(list.Index.kz)+1;
        v.aver = unique(list.Index.aver)+1;
        raw_params = list;
    otherwise
        v = -1;
end
end

%%read_cpx_header
function header = read_cpx_header(file, output)

%--------------------------------------------------------------------------
% header = read_cpx_header(file)
%
% read_cpx_header: Function for reading the header of a cpx file
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               output      can be 'yes' or 'no'. Default is 'yes'
%                           Specifies if the information about the cpx file is written
%                           to the command line
%
% Output:       header      gives out the header of the cpx file as an
%                           array. The structure of this array is the
%                           following:
%                           Column: 1:Stack, 2:Slice, 3:Coil, 4:Heart phase,
%                                   5:Echo, 6:Dynamics, 7:Segments,
%                                   8:data offset, 9:Scaling factor1,
%                                   10:Scaling factor2, 11:Compression,
%                                   12:Flip, 13: Scaling factor1, 14:
%                                   scaling factor2, 15: Mix, 16: Prep Dir,
%                                   17: Sequence Number 18: Segment2,
%                                   19: Syncho Nr.
%
%------------------------------------------------------------------------
header = [];

if nargin == 1
    output = 'yes';
end

% Calculate the size of the cpx-file
fid = fopen(file);
fseek(fid, 0, 'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');

% Read in the first header for
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 111,'long');

res_x = h1(11);
res_y = h1(12);
compression = h1(14);
if ~h2(26)
    offset = h1(10);
else
    offset = h2(26);
end

matrix_data_blocks = h1(13);

% Calculates the number of images in the cpx-file
image_exist = 1; i=0;
while image_exist
    %     header_offset = (res_x * res_y * 8 /compression + offset)*i;
    header_offset = (matrix_data_blocks * 512 + offset)*i;
    fseek(fid, header_offset, 'bof');
    h1 = fread(fid, 15, 'long');
    image_exist = h1(9);
    i = i+1;
end
images = i-1;

% Defines the header:
% header Columns : 1:Stack, 2:Slice, 3:Coil, 4:Heart phase, 5:Echo, 6:Dynamics,
%                  7:Segments, 8:data offset, 9:Resolution x, 10:Resolution y,
%                  11: Compression, 12: Flip, 13:Scaling factor1, 14:Scaling factor2
%                  15: Mix, 16: Prep Dir, 17: Sequence Nr.
%                  18: Segment2, 19: Syncho Number
header = zeros(images, 19);

% Runs over all images in the file and writes out its header
for i = 0: images-1
    header_offset = (matrix_data_blocks * 512 + offset)*i;
    fseek(fid, header_offset, 'bof');
    h1 = fread(fid, 15, 'long');
    factor = fread(fid,2,'float');
    h2 = fread(fid, 111,'long');
    header(i+1,1) = h1(2);                  % Stack
    header(i+1,2) = h1(3);                  % Slice
    header(i+1,3) = h2(2);                  % Coil
    header(i+1,4) = h1(6);                  % Heart phase
    header(i+1,5) = h1(5);                  % Echo
    header(i+1,6) = h1(7);                  % Dynamics
    header(i+1,7) = h1(8);                  % Segments
    if ~h2(26)
        header(i+1,8) = h1(10);                 % Data offset
    else
        header(i+1,8) = h2(26);                 % Data offset
    end
    header(i+1,9) = h1(11);                 % Resolution x
    header(i+1,10) = h1(12);                % Resolution y
    header(i+1,11) = h1(14);                % Compression
    header(i+1,12) = h2(111);               % Flip
    header(i+1,13) = factor(1);             % Scaling factor 1
    header(i+1,14) = factor(2);             % Scaling factor 2
    header(i+1,15) = h1(1);                 % mix
    header(i+1,16) = h1(4);                 % Prep Dir
    header(i+1,17) = h1(15);                % Sequence Number
    header(i+1,18) = h2(1);                 % Segment2
    header(i+1,19) = h2(3);                 % Syncho number
    
    if h1(9) == 0
        'Header Problem!! Too many images calculated'
        break
    end
end

% Reads in the last header and checks the parameter "Complex Matrix
% Existence" if it hasn't the value 0, the file is corrupt
last_header_offset = (matrix_data_blocks * 512 + offset)*images;
fseek(fid, last_header_offset, 'bof');
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 10,'long');
if h1(9) ~= 0
    'Header Problem'
    return
end

% Prints the parameters on the screen
if strcmp(output,'yes')
    s1=sprintf('\nResolution in x-direction: %d \nResolution in y-direction: %d \nNumber of stacks: %d \nNumber of slices: %d \nNumber of coils: %d \nNumber of heart phases: %d \nNumber of echos: %d \nNumber of dynamics: %d \nNumber of segments: %d \nNumber of segments2: %d',header(1,9),header(1,10),max(header(:,1))+1,max(header(:,2))+1,max(header(:,3))+1, max(header(:,4))+1,max(header(:,5))+1,max(header(:,6))+1,max(header(:,7))+1,max(header(:,18))+1);
    disp(s1);
end
end

%%create offset table
function offset_table_cpx=create_offset_table(header)

%--------------------------------------------------------------------------
% offset_table_cpx=create_offset_table(header)
%
% offset_table_cpx: Creates an array with all image offsets in the cpx-file
%                   The offset of a certain image can be obtained by:
%                   offset_table_cpx(stack, slice, coil, heart_phase, echo, dynamic, segment, segment2)
%
% Input:		header              The header of the cpx-file. Can be obtained
%                                   with the function "read_cpx_header"
%
% Output:       offset_table_cpx    An array with the image offsets in the
%                                   cpx-file
%
%
%------------------------------------------------------------------------


[rows, columns] = size(header);

for i = 1: rows
    if header(i,8) == 0;
        offset_table_cpx(header(i,1)+1, header(i,2)+1, header(i,3)+1, header(i,4)+1, header(i,5)+1, header(i,6)+1, header(i,7)+1, header(i,18)+1) = -100;
    else
        offset_table_cpx(header(i,1)+1, header(i,2)+1, header(i,3)+1, header(i,4)+1, header(i,5)+1, header(i,6)+1, header(i,7)+1, header(i,18)+1) = header(i,8);
    end
end
offset_table_cpx(find(offset_table_cpx==0)) = -1;
offset_table_cpx(find(offset_table_cpx==-100)) = 0;
end

%%read cpx image
function data = read_cpx_image(file, offset, border, flip_img)

%--------------------------------------------------------------------------
% data = read_cpx_image(file, offset, border, flip_img)
%
% read_cpx_image : Function for reading one image, defined through the
%                  input parameters
%
% Input:		file		cpx-file must be given in ''. e.g. 'file.cpx'
%               offset      offset to the image from header (element (.,8))
%               border      0 returns the original sized image
%                           1 makes a square image
%               flip_img    0 no flip of the image
%                           1 flips the image (like the scanner does)
%
%               Notice! All numbers starts from 1. E.g.the first slice is
%               numbered with 1
%
% Output:       data        The requested image in a 2 dimensional array
%
%
% Notice! For scans with body coil and surface coil in different stacks,
% e.g. SENSE reference scan, the body coil is numbered as coil no. 1 in
% stack no. 2, the surface coils are numbered as coils 2,3,4... in stack no. 1!
%------------------------------------------------------------------------


%Reads the header of the requested image
fid = fopen(file);

fseek(fid, offset-512,'bof');
h1 = fread(fid, 15, 'long');
factor = fread(fid,2,'float');
h2 = fread(fid, 10,'long');

res_x = h1(11);
res_y = h1(12);
compression = h1(14);

%Reads the requested image
fseek(fid, offset,'bof');
switch (compression)
    case 1
        data = zeros(res_x*res_y*2,1,'single');
        data=fread(fid, res_x*res_y*2, 'float');
    case 2
        data = zeros(res_x*res_y*2,1,'single');
        data(:)=fread(fid, res_x*res_y*2,'short');
        data=factor(2)+factor(1).*data;
    case 4
        data = zeros(res_x*res_y*2,1,'single');
        data=fread(fid, res_x*res_y*2, 'int8');
        data=factor(2)+factor(1).*data;
end
data = complex(data(1:2:end),data(2:2:end));
data = reshape(data,res_x,res_y);

%Adds the border if requested
if border & (res_x ~= res_y)
    res = max([res_x, res_y]);
    data_temp = zeros(res, res);
    if res_x > res_y
        data_temp(:,floor((res - res_y)/2): res - ceil((res - res_y)/2+0.1)) = data;
    else
        data_temp(floor((res - res_x)/2): res - ceil((res - res_x)/2+0.1),:) = data;
    end
    data = data_temp;
    clear data_temp;
end

%Flips the image if requested
if flip_img
    s = size(data);
    data = data(end:-1:1,:);
    %     data=data';
    %     data = data(:,s(1)+1-(1:s(1)));
    %     data = data(s(2)+1-(1:s(2)),:);
end

fclose(fid);
end

%%load noise from raw file
function [data,info] = read_noise(filename,varargin)

% Start execution time clock and initialize DATA and INFO to empty arrays
tic;
data=[];
info=[];

% Initialize INFO structure
% Serves to fix the display order
info.filename = [];
info.loadopts = [];
info.dims = [];
info.labels = [];
info.labels_row_index_array = [];
info.label_fieldnames = [];
info.idx = [];
info.fseek_offsets = [];
info.nLabels = [];
info.nLoadedLabels = [];
info.nDataLabels = [];
info.nNormalDataLabels = [];
info.datasize = [];


% Parse the filename.
% It may be the LAB filename, RAW filename or just the filename prefix
% Instead of REGEXP, use REGEXPI which igores case
toks = regexpi(filename,'^(.*?)(\.lab|\.raw)?$','tokens');
prefix = toks{1}{1};
labname = sprintf('%s.lab',prefix);
rawname = sprintf('%s.raw',prefix);
info.filename = filename;
% Open LAB file and read all hexadecimal labels
labfid = fopen(labname,'r');
if labfid == -1
    error('Cannot open %s for reading', labname);
end

% Read all hexadecimal labels
unparsed_labels = fread (labfid, [16 Inf], 'uint32=>uint32');
info.nLabels = size(unparsed_labels,2);
fclose(labfid);

% Parse hexadecimal labels
% Inspired by Holger Eggers' readRaw.m.  Thanks Holger! 
% See arsrcglo1.h for more details.
info.labels.DataSize.vals         = unparsed_labels(1,:);
info.labels.LeadingDummies.vals   = bitshift (bitand(unparsed_labels(2,:), (2^16-1)),  -0);
info.labels.TrailingDummies.vals  = bitshift (bitand(unparsed_labels(2,:), (2^32-1)), -16);
info.labels.SrcCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^16-1)),  -0);
info.labels.DstCode.vals          = bitshift (bitand(unparsed_labels(3,:), (2^32-1)), -16);
info.labels.SeqNum.vals           = bitshift (bitand(unparsed_labels(4,:), (2^16-1)),  -0);
info.labels.LabelType.vals        = bitshift (bitand(unparsed_labels(4,:), (2^32-1)), -16);
info.labels.ControlType.vals      = bitshift( bitand(unparsed_labels(5,:),  (2^8-1)),  -0);
info.labels.MonitoringFlag.vals   = bitshift( bitand(unparsed_labels(5,:), (2^16-1)),  -8);
info.labels.MeasurementPhase.vals = bitshift( bitand(unparsed_labels(5,:), (2^24-1)), -16);
info.labels.MeasurementSign.vals  = bitshift( bitand(unparsed_labels(5,:), (2^32-1)), -24);
info.labels.GainSetting.vals      = bitshift( bitand(unparsed_labels(6,:),  (2^8-1)),  -0);
info.labels.Spare1.vals           = bitshift( bitand(unparsed_labels(6,:), (2^16-1)),  -8);
info.labels.Spare2.vals           = bitshift (bitand(unparsed_labels(6,:), (2^32-1)), -16);
info.labels.ProgressCnt.vals      = bitshift (bitand(unparsed_labels(7,:), (2^16-1)),  -0);
info.labels.Mix.vals              = bitshift (bitand(unparsed_labels(7,:), (2^32-1)), -16);
info.labels.Dynamic.vals          = bitshift (bitand(unparsed_labels(8,:), (2^16-1)),  -0);
info.labels.CardiacPhase.vals     = bitshift (bitand(unparsed_labels(8,:), (2^32-1)), -16);
info.labels.Echo.vals             = bitshift (bitand(unparsed_labels(9,:), (2^16-1)),  -0);
info.labels.Location.vals         = bitshift (bitand(unparsed_labels(9,:), (2^32-1)), -16);
info.labels.Row.vals              = bitshift (bitand(unparsed_labels(10,:), (2^16-1)),  -0);
info.labels.ExtraAtrr.vals        = bitshift (bitand(unparsed_labels(10,:), (2^32-1)), -16);
info.labels.Measurement.vals      = bitshift (bitand(unparsed_labels(11,:), (2^16-1)),  -0);
info.labels.E1.vals               = bitshift (bitand(unparsed_labels(11,:), (2^32-1)), -16);
info.labels.E2.vals               = bitshift (bitand(unparsed_labels(12,:), (2^16-1)),  -0);
info.labels.E3.vals               = bitshift (bitand(unparsed_labels(12,:), (2^32-1)), -16);
info.labels.RfEcho.vals           = bitshift (bitand(unparsed_labels(13,:), (2^16-1)),  -0);
info.labels.GradEcho.vals         = bitshift (bitand(unparsed_labels(13,:), (2^32-1)), -16);
info.labels.EncTime.vals          = bitshift (bitand(unparsed_labels(14,:), (2^16-1)),  -0);
info.labels.RandomPhase.vals      = bitshift (bitand(unparsed_labels(14,:), (2^32-1)), -16);
info.labels.RRInterval.vals       = bitshift (bitand(unparsed_labels(15,:), (2^16-1)),  -0);
info.labels.RTopOffset.vals       = bitshift (bitand(unparsed_labels(15,:), (2^32-1)), -16);
info.labels.ChannelsActive.vals   = unparsed_labels(16,:);

clear unparsed_labels;

% Find unique values of each label field
info.label_fieldnames = fieldnames(info.labels);
for k=1:length(info.label_fieldnames),
    info.labels.(info.label_fieldnames{k}).uniq = unique( info.labels.(info.label_fieldnames{k}).vals ); 
end

% Calculate fseek offsets
info.fseek_offsets = zeros(info.nLabels,1);
info.fseek_offsets(1)=512; % add mysterious 512 byte offset to begin reading file
for k=2:info.nLabels,
    info.fseek_offsets(k) = info.fseek_offsets(k-1)+ info.labels.DataSize.vals(k-1) - info.labels.TrailingDummies.vals(k-1)  - info.labels.LeadingDummies.vals(k-1);
end
info.idx.no_data = find(info.labels.DataSize.vals==0);
info.fseek_offsets(info.idx.no_data) = -1;

% Find indices of different label control types
% See arsrcglo1.h for more details.
standard_labels = info.labels.LabelType.vals==32513;
info.idx.NORMAL_DATA         = find(info.labels.ControlType.vals== 0 & standard_labels);
info.idx.DC_OFFSET_DATA      = find(info.labels.ControlType.vals== 1 & standard_labels);
info.idx.JUNK_DATA           = find(info.labels.ControlType.vals== 2 & standard_labels);
info.idx.ECHO_PHASE_DATA     = find(info.labels.ControlType.vals== 3 & standard_labels);
info.idx.NO_DATA             = find(info.labels.ControlType.vals== 4 & standard_labels);
info.idx.NEXT_PHASE          = find(info.labels.ControlType.vals== 5 & standard_labels);
info.idx.SUSPEND             = find(info.labels.ControlType.vals== 6 & standard_labels);
info.idx.RESUME              = find(info.labels.ControlType.vals== 7 & standard_labels);
info.idx.TOTAL_END           = find(info.labels.ControlType.vals== 8 & standard_labels);
info.idx.INVALIDATION        = find(info.labels.ControlType.vals== 9 & standard_labels);
info.idx.TYPE_NR_END         = find(info.labels.ControlType.vals==10 & standard_labels);
info.idx.VALIDATION          = find(info.labels.ControlType.vals==11 & standard_labels);
info.idx.NO_OPERATION        = find(info.labels.ControlType.vals==12 & standard_labels);
info.idx.DYN_SCAN_INFO       = find(info.labels.ControlType.vals==13 & standard_labels);
info.idx.SELECTIVE_END       = find(info.labels.ControlType.vals==14 & standard_labels);
info.idx.FRC_CH_DATA         = find(info.labels.ControlType.vals==15 & standard_labels);
info.idx.FRC_NOISE_DATA      = find(info.labels.ControlType.vals==16 & standard_labels);
info.idx.REFERENCE_DATA      = find(info.labels.ControlType.vals==17 & standard_labels);
info.idx.DC_FIXED_DATA       = find(info.labels.ControlType.vals==18 & standard_labels);
info.idx.DNAVIGATOR_DATA     = find(info.labels.ControlType.vals==19 & standard_labels);
info.idx.FLUSH               = find(info.labels.ControlType.vals==20 & standard_labels);
info.idx.RECON_END           = find(info.labels.ControlType.vals==21 & standard_labels);
info.idx.IMAGE_STATUS        = find(info.labels.ControlType.vals==22 & standard_labels);
info.idx.TRACKING            = find(info.labels.ControlType.vals==23 & standard_labels);
info.idx.FLUOROSCOPY_TOGGLE  = find(info.labels.ControlType.vals==24 & standard_labels);
info.idx.REJECTED_DATA       = find(info.labels.ControlType.vals==25 & standard_labels);
info.idx.UNKNOWN27           = find(info.labels.ControlType.vals==27 & standard_labels);
info.idx.UNKNOWN28           = find(info.labels.ControlType.vals==28 & standard_labels);

% Calculate number of standard, normal data labels
info.nNormalDataLabels = length(info.idx.NORMAL_DATA);

% Dimension names
dimnames = {'coil','kx','ky','kz','E3','loc','ec','dyn','ph','row','mix','avg'};
dimfields = {'N/A','N/A','E1','E2','E3','Location','Echo','Dynamic','CardiacPhase','Row','Mix','Measurement'};

% Initialize dimension data to zero
info.dims.nCoils         = 0;
info.dims.nKx            = 0;
info.dims.nKy            = 0;
info.dims.nKz            = 0;
info.dims.nE3            = 0;
info.dims.nLocations     = 0;
info.dims.nEchoes        = 0;
info.dims.nDynamics      = 0;
info.dims.nCardiacPhases = 0;
info.dims.nRows          = 0;
info.dims.nMixes         = 0;
info.dims.nMeasurements  = 0;

% Calculate max number of active coils
maxChannelsActiveMask = 0;
for k=1:length(info.labels.ChannelsActive.uniq),
    maxChannelsActiveMask = bitor(maxChannelsActiveMask,info.labels.ChannelsActive.uniq(k));
end
while maxChannelsActiveMask > 0
    if bitand(maxChannelsActiveMask, 1),
        info.dims.nCoils = info.dims.nCoils + 1;
    end
    maxChannelsActiveMask = bitshift (maxChannelsActiveMask, -1);
end

% Calculate dimensions of normal data
info.dims.nKx            = max(info.labels.DataSize.vals(info.idx.NORMAL_DATA)) / info.dims.nCoils / 2 / 2;
info.dims.nKy            = length(unique(info.labels.E1.vals(info.idx.NORMAL_DATA)));
info.dims.nKz            = length(unique(info.labels.E2.vals(info.idx.NORMAL_DATA)));
info.dims.nE3            = length(unique(info.labels.E3.vals(info.idx.NORMAL_DATA)));
info.dims.nLocations     = length(unique(info.labels.Location.vals(info.idx.NORMAL_DATA)));
info.dims.nEchoes        = length(unique(info.labels.Echo.vals(info.idx.NORMAL_DATA)));
info.dims.nDynamics      = length(unique(info.labels.Dynamic.vals(info.idx.NORMAL_DATA)));
info.dims.nCardiacPhases = length(unique(info.labels.CardiacPhase.vals(info.idx.NORMAL_DATA)));
info.dims.nRows          = length(unique(info.labels.Row.vals(info.idx.NORMAL_DATA)));
info.dims.nMixes         = length(unique(info.labels.Mix.vals(info.idx.NORMAL_DATA)));
info.dims.nMeasurements  = length(unique(info.labels.Measurement.vals(info.idx.NORMAL_DATA)));

% With known possible dimension names, the load options can now be parsed
p = inputParser;
p.StructExpand = true;
p.CaseSensitive = true;
p.KeepUnmatched = false; % throw an error for unmatched inputs
p.addRequired('filename', @ischar);
for k=1:length(dimnames)
    p.addParameter(dimnames{k}, [], @isnumeric);
end
p.addParameter('verbose', false, @islogical);
p.addParameter('savememory', true, @islogical);
p.parse(filename, varargin{:});

% Return loadopts structure inside INFO structure
% remove filename field - it is passed as the first required argument
info.loadopts = rmfield(p.Results,'filename');

% Find the unique set of values for each dimension name
info.dims.coil = 1:info.dims.nCoils;
info.dims.kx   = 1:info.dims.nKx;
for k=3:length(dimnames) % skip coil and kx
    info.dims.(dimnames{k}) = unique(info.labels.(dimfields{k}).vals(info.idx.NORMAL_DATA));
end

% Find intersection of available dimensions with LOADOPTS dimensions
for k=1:length(dimnames)
    if ~isempty(info.loadopts.(dimnames{k}))
        info.dims.(dimnames{k}) = intersect_a_with_b(info.loadopts.(dimnames{k}),info.dims.(dimnames{k}));
    end
end

% Calculate data size
datasize = []; 
for k=1:length(dimnames)
    datasize = [datasize length(info.dims.(dimnames{k}))];
end
info.datasize = datasize;

% throw error if any dimension size is zero
if any(info.datasize==0)
    zero_length_str = sprintf(' ''%s'' ', dimnames{info.datasize==0});
    error('size of selected data to load has zero length along dimension(s): %s', zero_length_str);
end

% Skip data loading if only one output argument is provided, return INFO
if nargout==1
    info.labels_row_index_array = 1:size(info.labels,1);
    data=info;
    return;
end

% Create array to hold label row numbers for loaded data
% skip the coil and kx dimensions
info.labels_row_index_array = zeros(datasize(3:end));

% Pre-allocate DATA array
if info.loadopts.savememory==true
    data = zeros(info.datasize,'single');
else
    data = zeros(info.datasize);
end

% Read RAW data for selected dimension ranges
fidraw = fopen(rawname,'r','ieee-le');
if fidraw == -1
    error('cannot open RAW file: %s', rawname);
end
info.nLoadedLabels=0;

raw_data_fread_size = double(info.dims.nCoils * info.dims.nKx * 2);
rawdata_2d = complex(zeros(info.dims.nCoils,info.dims.nKx),zeros(info.dims.nCoils,info.dims.nKx));

% Read FRC noise data
if(~isempty(info.idx.FRC_NOISE_DATA))
    
    for n=1:length(info.idx.FRC_NOISE_DATA)
        
        frc_noise_idx = info.idx.FRC_NOISE_DATA(n);
        
        num = double(info.labels.ChannelsActive.vals(frc_noise_idx));
        [f,e]=log2(num);
        ncoils = sum(rem(floor(num*pow2(1-e:0)),2));
        
        frc_noise_samples_per_coil = info.labels.DataSize.vals(frc_noise_idx) / 2 / 2 / ncoils;
    byte_offset = info.fseek_offsets(frc_noise_idx);
    status = fseek(fidraw, byte_offset, 'bof');
    rawdata_1d = double(fread(fidraw, double(info.labels.DataSize.vals(frc_noise_idx)/2) , 'int16'));
    info.FRC_NOISE_DATA(1:ncoils,:,n) = permute(reshape(complex(rawdata_1d(1:2:end), rawdata_1d(2:2:end)), frc_noise_samples_per_coil, ncoils),[2 1]);
    end
end
fclose(fidraw);

% Calculate total raw data blobs
size_data = size(data);
max_img_dims = size_data(3:end);
info.nDataLabels = prod(max_img_dims);

% If VERBOSE, display execution information
if info.loadopts.verbose==true
    disp( sprintf('Loaded %d of %d available normal data labels', info.nLoadedLabels, info.nNormalDataLabels) );
    tmpstr = '';
    for k=1:length(dimnames)
        tmpstr = sprintf('%s, # %s: %d', tmpstr, dimnames{k}, length(info.dims.(dimnames{k})) );
    end
    disp( sprintf('Data contains %d raw labels - %s', info.nDataLabels, tmpstr(3:end)) );
    fprintf('Total execution time = %.3f seconds\n', toc) ;
end

% Find intersection of vector a with vector b without sorting 
function c = intersect_a_with_b(a,b)
c = a;
% work backwards in order to use [] assignment
for k=length(a):-1:1
    if isempty(find(a(k)==b,1))
        c(k)=[]; 
    end
end

% force c to be a row vector for easier display
c = c(:).';


end
end

%%water removal
function FIDsvd = waterremovalSVD(FID, sw, nHSVD, HSVDLow, HSVDHigh, plottingon)
%***********************************************************
% waterSVD.m
%
% Residual water resonance removal by SVD.
%
% 1. SVD of spectral region [ZxLow ... ZxHigh]
% 2. Least-squares fitting to model the total spectrum
% 3. Extraction of resonances in the region [WxLow ... WxHigh]
% 4. Subtraction of water from original FID
%
% By Robin A. de Graaf
% MRRC, Yale University
% Original : March, 1999
% Modified : September, 2008
%***********************************************************
if plottingon
    disp('starting SVD water removal.....')
end

DataHandling = 3;

zff = 1;
% global amp frq decay phs;

% clear i x H U V S Utr Vtr Str;

FIDorig = FID;
clear FID;

nporig = length(FIDorig(:));
npsvd = nporig/2; % VB take half the points for speed

FID = FIDorig(1:npsvd);

% f= length(FID);			% Number of complex datapoints

FIDsvd = FID(1:np);

npnonzero = length(FIDsvd);

tacq = npnonzero/(1000*sw);                % Acquisition time (in s)
dt = tacq/npnonzero;                       % Dwell-time (in s)
time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
time = reshape(time,npnonzero,1);

Lmax = round(0.4*npnonzero);               % Dimension 1 for LxM SVD matrix
Mmax = npnonzero+1-Lmax;                   % Dimension 2 for LxM SVD matrix

if plottingon
    disp(' ');
    disp('Water removal in progress ... ');
end

%*****************
% Allocate memory
%*****************
U = zeros(Lmax,Lmax);
V = zeros(Mmax,Mmax);
S = zeros(Lmax,Mmax);

H = zeros(Lmax,Mmax);

%**********************************************
% Create Hankel matrix from original FID data
%**********************************************
for L = 1:1:Lmax;
    M = 1:1:Mmax;
    H(L,M) = FIDsvd(L+M-1);
end;

%**********************************************
% Perform SVD on Hankel matrix
%**********************************************
if plottingon
    tic;
    disp('Step 1 : SVD of Hankel matrix in progress ...');
end

[U,S,V] = svd(H);

if plottingon
    tt = toc;
    dd = ['... done in ' num2str(tt,3) ' s.'];
    disp(dd);
end

Utr = zeros(Lmax,nHSVD);
Vtr = zeros(nHSVD,Mmax);
Str = zeros(nHSVD,nHSVD);

Uup = zeros(Lmax-1,nHSVD);
Udown = zeros(Lmax-1,nHSVD);
Udownc = zeros(nHSVD,Lmax-1);

Z = zeros(nHSVD,nHSVD);
W = zeros(nHSVD,nHSVD);
D = zeros(nHSVD,nHSVD);

if plottingon
    disp('Step 2 : Calculation of lineshape parameters in progress ...');
end

tic; 

%**********************************************
% Calculate truncated SVD matrix
%**********************************************
Utr = U(:,1:1:nHSVD);
Str = S(:,1:1:nHSVD);
Vtr = V(:,1:1:nHSVD);

for kk1 = 2:Lmax;
   for kk2 = 1:nHSVD;
      Uup(kk1-1,kk2) = Utr(kk1,kk2);
      Udown(kk1-1,kk2) = Utr(kk1-1,kk2);
   end;
end;

Z = pinv(Udown)*Uup;

q = eig(Z);
q = log(q);

%******************************************************
% Determination of frequencies and T2 constants from D
%******************************************************
clear frq decay;

% Frequency (in Hz)
frq = imag(q)/(2*pi*dt);
% Time constant (in s)
decay = real(q)/dt;

clear ampcomplex amp phs basis;

switch DataHandling
    case 1
        time = 0:dt:(npnonzero-1)*dt;              % Time base of FID
        time = reshape(time,npnonzero,1);

        % Calculate basis functions
        for kk1 = 1:1:nHSVD;
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end;
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDsvd;
    case 2
        time = 0:dt:(npsvd-1)*dt;              % Time base of FID
        time = reshape(time,npsvd,1);

        % Calculate basis functions
        for kk1 = 1:1:nHSVD;
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end;
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDsvd;
    case 3
        time = 0:dt:(nporig-1)*dt;              % Time base of FID
        time = reshape(time,nporig,1);

        % Calculate basis functions
        for kk1 = 1:1:nHSVD;
            basis(:,kk1) = exp((decay(kk1)+2*pi*1i*frq(kk1))*time);
        end;
        
        % Amplitude estimates
        ampcomplex = pinv(basis)*FIDorig;
end;

amp = abs(ampcomplex);
phs = atan2(imag(ampcomplex),real(ampcomplex));

if plottingon
    tt = toc;
    dd1 = [' ... done in ' num2str(tt,3) ' s.'];
    disp(dd1);

    disp('Step 3 : Spectral reconstruction in progress ...');

    tic; 
end

%***********************************************************************
% Reconstruct signal in the frequency range [HSVDLow, HSVDHigh]
%***********************************************************************
waterpos = find(((frq > -1000*HSVDHigh) & (frq < -1000*HSVDLow)));

nwater = length(waterpos);

switch DataHandling
    case 1
        FIDw = zeros(npnonzero,1);
        FIDw0 = zeros(npnonzero,1);
    case 2
        FIDw = zeros(npsvd,1);
        FIDw0 = zeros(npsvd,1);
    case 3
        FIDw = zeros(nporig,1);
        FIDw0 = zeros(nporig,1);
end;

meanfrq = mean(frq(waterpos(1:nwater)));

for kk1 = 1:1:nwater;
   FIDcomponent = amp(waterpos(kk1)).*exp(2*pi*1i*frq(waterpos(kk1))*time).*exp(time*decay(waterpos(kk1))).*exp(1i*phs(waterpos(kk1)));
   FIDcomponent0 = amp(waterpos(kk1)).*exp(2*pi*1i*(frq(waterpos(kk1))-meanfrq)*time).*exp(time*decay(waterpos(kk1))).*exp(1i*phs(waterpos(kk1)));
   FIDw = FIDw + FIDcomponent;
   FIDw0 = FIDw0 + FIDcomponent0;
   clear FIDcomponent FIDcomponent0;
end;

switch DataHandling
    case 1
        if (npsvd ~= npnonzero)
            FIDw(npsvd+1:npsvd) = 0;
            FIDw0(npsvd+1:npsvd) = 0;
            time = 0:dt:(npsvd-1)*dt;              % Time base of FID
            tacq = (npsvd-1)*dt;
        end;
    case 2
        if (npsvd ~= npnonzero)
            FIDw(npsvd+1:nporig) = 0;
            FIDw0(npsvd+1:nporig) = 0;
            time = 0:dt:(nporig-1)*dt;              % Time base of FID
            tacq = (nporig-1)*dt;
        else
            FIDw(npsvd+1:nporig) = 0.0;
            FIDw0(npsvd+1:nporig) = 0.0;
            time = 0:dt:(nporig-1)*dt;              % Time base of FID
            tacq = (nporig-1)*dt;
        end;
        FID = FIDorig;
    case 3
        FID = FIDorig;
end;

% Fourier transformation
spec = fftshift(fft(FIDorig,zff*nporig));
spec = reshape(spec,length(spec),1);
specw = fftshift(fft(FIDw,zff*nporig));
specw = reshape(specw,length(specw),1);

% Phase correction
% spec = SVDPhaseCorrection(spec, ph0, ph1, sw);
specA = real(spec);
% specw = SVDPhaseCorrection(specw, ph0, ph1, sw);
specAw = real(specw);

if plottingon
    tt = toc;
    dd = ['... done, using ' num2str(nwater) ' water components, in ' num2str(tt,3) ' s.'];
    disp(dd);
    disp('Water removal completed.');
end

%*********************************************************************
% Display original and fitted FID and spectra and the differences
%*********************************************************************
if plottingon
    hh = figure(4);
    set(hh,'position',[200 50 800 600])

    freq = 0.5*sw:(-sw/(zff*nporig-1)):-0.5*sw;

    subplot(2,2,1), plot(time,real(FID),'b',time,real(FIDw),'r');
    axis([0 tacq 1.1*min(real(FID)) 1.1*max(real(FID))])
    title('Original (blue)/fitted (red) FID');
    subplot(2,2,3), plot(time,real(FID-FIDw));
    axis([0 tacq 1.1*min(real(FID-FIDw)) 1.1*max(real(FID-FIDw))])
    title('Difference FID');

    subplot(2,2,2), plot(freq,specA,'b',freq,specAw,'r');
%     axis([ZxLow ZxHigh ZyLow ZyHigh])
    title('Original (blue)/fitted (red) spectrum');
    subplot(2,2,4), plot(freq,(specA-specAw));
%     axis([ZxLow ZxHigh ZyLow ZyHigh])
    title('Difference spectrum');
end

FIDsvd = FID - FIDw;

end


% write sdat spar
function write_sdat_spar(data, fname_p, fname_out, fname_out_spar, averages, np)

data = data(:);
data1 = zeros(size(data,1)*2,1);
data1(1:2:end) = real(data);
data1(2:2:end) = imag(data);

fid2 = fopen(fname_out,'w','ieee-le');
if fid2 ~= -1
    disp(['storing file: ' fname_out])
    fwriteVAXG(fid2,data1,'float32');
    fclose(fid2);
    
    fid4 = fopen(fname_p,'r');
    fid5 = fopen(fname_out_spar,'w');
    
    while ~feof(fid4);
        a = fgetl(fid4);
        if ~isempty(strfind(a,'samples : '))
            fprintf(fid5,'samples : %d\r\n',np);
        elseif ~isempty(strfind(a,'rows : '))
            fprintf(fid5,'rows : %d\r\n',averages);
        elseif ~isempty(strfind(a,'averages : '))
            fprintf(fid5,'averages : %d\r\n',1);
        elseif ~isempty(strfind(a,'spec_num_row : '))
            fprintf(fid5,'spec_num_row : %d\r\n',averages);
        elseif ~isempty(strfind(a,'dim2_pnts : '))
            fprintf(fid5,'dim2_pnts : %d\r\n',averages);
        else
            fprintf(fid5,'%s\r\n',a);
        end
        
    end
    fclose(fid4);
    fclose(fid5);
%     copyfile(fname_p,fname_out_spar,'f');
else
    disp(['couldnt open ' fname_out])
end

end

function count = fwriteVAXG(fid, A, precision)
% FWRITEVAXG(FID, A , PRECISION) writes the elements of 'A' to the 
% specified file in VAXG format, translating MATLAB ® values to the specified
% precision. COUNT is the number of elements successfully written.  
% 
% FID is an integer file identifier obtained from FOPEN. FWRITEVAXG requires
% FOPEN to open the output file in IEEE little-endian machineformat.
%
% PRECISION controls the form and size of result. See FREAD for supported
% formats.
%
% Usage:
%   A = rand(3,3); 
% 	fid = fopen('myFile', 'w', 'ieee-le');
%	count = fwriteVAXG(fid, A, 'double');
%   fclose(fid);
%
% The function is intended to be called by the user.
%
% ® 2009 The MathWorks, Inc. MATLAB and Simulink are registered trademarks
% of The MathWorks, Inc. See www.mathworks.com/trademarks for a list of 
% additional trademarks. Other product or brand names may be trademarks or 
% registered trademarks of their respective holders.

% Check for proper number of input arguments
if nargin < 2
    error('Not enough input arguments.')
end

% Check that the file has been open in ieee-le machineformat
[~, ~, machineformat] = fopen(fid);
if ~strcmp(machineformat, 'ieee-le')
    error('Use FOPEN with ieee-le precision');
end

switch precision

    case {'float32', 'single'}
      rawUINT32 = VAXF_to_uint32le(A);
      count = fwrite(fid, rawUINT32, 'uint32');

    case {'float64', 'double'}
      rawUINT32 = VAXG_to_uint64le(A);
      count = fwrite(fid, rawUINT32, 'uint32');
      count = count/2;%2 UINT32 pieces for each double precision number

    case {'float'}
      if intmax == 2147483647 %32bit OS float is 32 bits
         rawUINT32 = VAXF_to_uint32le(A);
         count = fwrite(fid, rawUINT32, 'uint32');
      else
         rawUINT32 = VAXG_to_uint64le(A);
         count = fwrite(fid, rawUINT32, 'uint32');
         count = count/2;%2 UINT32 pieces for each double precision number
      end

    otherwise

      count = fwrite(fid, A, precision, 'vaxg');

end

end

function [ uint32le] = VAXF_to_uint32le(floatVAXF)
%VAXF_TO_UINT32LE Converts from VAXF (single precision) to IEEE-LE (UINT32)


A = 2;      % VAX specific
B = 128;    % VAX specific
C = 0.5;    % VAX specific
D = log(2); % VAX specific

% Determine the sign bit. If -ve transform to positive.
S = zeros(size(floatVAXF));
if any(floatVAXF(:) < 0)
    indices = find(floatVAXF<0);
    floatVAXF(indices) = (-1) .* floatVAXF(indices);
    S = zeros(size(floatVAXF));
    S(indices) = 1;
end

% Decompose the floating point number to SEF (Sign, Exp, Fract)
E = floor((log(floatVAXF)./ D) + 1 + B); 
F = ((floatVAXF ./ A.^(double(E)-B))) - C; 
% Convert floating point fraction to unsigned integer
F = floor(F * 16777216);   %VAX Specific 16777216=2^24

% Shift the bits of S, E and F
S = bitshift(bitshift(uint32(S),0), 31);
E = bitshift(bitshift(uint32(E),0), 24);
F = bitshift(bitshift(uint32(F),0), 9);

% Combine the S, E and F into the unsigned integer value
vaxInt = bitor(bitor(S,bitshift(E, -1)),bitshift(F,-9));

% Swap WORD1 and WORD2
% VAX      <-----WORD1-----><-----WORD2----->
% IEEE-LE  <-----WORD2-----><-----WORD1----->

word1 = bitshift(vaxInt,16);
word2 = bitshift(vaxInt,-16);

uint32le = bitor(word1,word2);

	
end
