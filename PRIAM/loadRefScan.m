function [ref_file,Ref,img_array,img_body,noise_array,noise_body] = loadRefScan(spec_path)
%% [ref_file,Ref,img_array,img_body,noise_array,noise_body] = loadRefScan(spec_path)
%   Loads a Philips coil reference scan (*.cpx format) when provided a
%   path. Assumes a subdirectory '/ref' containing this data, relative to
%   the MRS data path.
%
%   Input:
%       spec_path   Folder containing the MRS data. loadRefScan assumes a
%                   subdirectory '/ref' containing the *.cpx data.
%
%   Output:
%       ref_file    cpx file used
%       Ref         struct containing dimensional information
%       img_array   4D stack of receiver coil images
%                   [number of coils, x, y, z]
%       img_body    4D stack of body coil images
%       noise_array 2D stack of receiver coil noise levels
%                   [number of coils, number of pixels estimated for noise]
%       noise_body  2D stack of body coil noise levels
%
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-03-15)
%       goeltzs1@jhmi.edu
%   
%   Credits:    
%       This code is based on an initial PRIAM reconstruction routine.
%       Dr. Vincent O. Boer (vincentob@drcmr.dk)
%       Danish Research Centre for Magnetic Resonance (Hvidovre Hospital)
%
%   History:
%       2018-03-15: First version of the code.
%


%% Put together the filenames for the coil reference scan

% Coil reference scan needs to be a .cpx file in a subfolder '/ref'
% relative to the MRS data.
disp('Determining coil reference scan...')
cpx_files = dir([spec_path filesep 'ref' filesep '*.cpx']);
% If more than one .cpx file is found, select one
if(length(cpx_files)>1)
    disp('More than one .cpx file found:');
    for ii = 1:length(cpx_files)
        fprintf('# %i --- %s\n',ii,cpx_files(ii).name);
    end
    result = input('Select the correct coil reference scan (input number): ');
    if isempty(result)
        return
    end
else
    result = 1;
end
ref_file = [spec_path filesep 'ref' filesep cpx_files(result).name];

%% Read coil reference imaging scan (*.cpx)

% Load the selected cpx file
fprintf('Loading coil reference scan...\n%s',ref_file);
img = read_cpx(ref_file);

% Permute the coil reference image array to be in the following order:
% 1 - number of coils
% 2 - image dimension x
% 3 - image dimension y
% 4 - image dimension z
% 5 - number of stacks:
%   1st stack contains receiver coil images
%   2nd stack contains body coil images (and is empty if less body coil
%   elements than receiver coil images)
img = permute(img,[5 1 2 4 3]);
Ref.ncoils = size(img,1);
Ref.dimx = size(img,2);
Ref.dimy = size(img,3);
Ref.dimz = size(img,4);

% Separate into two separate arrays
img_array = double(img(:,:,:,:,1));
img_body = double(img(:,:,:,:,2));

% Normalize within each array
img_array = img_array/max(abs(img_array(:)));
img_body = img_body/max(abs(img_body(:)));

% If there is a *.raw file for the reference image, extract the
% noise from the header.
if exist([ref_file(1:end-4) '.raw'],'file')
    disp('Found *.raw file for the coil reference scan. Extracting noise from the *.raw header.');
    [~, info] = read_noise([ref_file(1:end-4) '.raw']);
    noise_array = info.FRC_NOISE_DATA(:,:,1);
    noise_body = info.FRC_NOISE_DATA(:,:,2);
else
    % If not, extract noise levels from the edge of the first slice.
    disp('No *.raw file for the coil reference scan found. Extracting noise from the edges of the image.');
    noise_array = squeeze(img_array(:,1:30,:,1));
    noise_array = reshape(noise_array,[Ref.ncoils size(noise_array,2)*size(noise_array,3)]);
    noise_body = squeeze(img_body(:,1:30,:,1));
    noise_body = reshape(noise_body,[Ref.ncoils size(noise_body,2)*size(noise_body,3)]);
end

save([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_sense_coil_img.mat'],'ref_file','Ref','img_array','noise_array','noise_body');
save([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_volume_coil_img.mat'],'ref_file','Ref','img_body','noise_array','noise_body');
disp('Loading coil reference scan finished!');