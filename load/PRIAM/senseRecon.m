function MRS_struct = senseRecon(MRS_struct)
%% MRS_struct = senseRecon(MRS_struct)
%   Reads Siemens TWIX files (*.dat) and removes participant information. New
%   de-identified TWIX files are then output, with filenames appended with
%   '_noID'. The original files are not overwritten.
%
%   Input:
%       TWIXDeIdentify, by itself, de-identifies all TWIX files found
%       within the current directory.
%
%   Output:
%       c = {'MRS_01.dat', 'MRS_02.dat'};
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

%% Setup paths and filenames, and load coil reference scan

% Compute coil sensitivities from reference imaging scan, or from water MRS
% scan?
% 0 = reference imaging scan
% 1 = water MRS scan
MRS_struct.p.SENSE.sens_from_ref_scan = 0; % 1 not implemented yet, GO 11/01/2016

% Set path and filenames
spec_path = pwd;
% spec_file = MRS_struct.metabfile;
% spec_file = [spec_path filesep spec_file];
% spec_file = 'raw_priam_trans_neg.data';
% spec_path = '/Users/Georg/Documents/WORK_BALTIMORE/UserRequests/Gannet/2018/Leuven_MEGAPRIAM/03/transneg';

% Generate new folder GannetRecon_output
if ~exist([spec_path filesep 'GannetRecon_output'],'dir')
    mkdir([spec_path filesep 'GannetRecon_output']);
end

% Ask for spatial separation of the voxels that was entered into the exam
% card (in mm, positive value)
MRS_struct.p.vox_sep = input('What is the separation of the voxels in mm?  ');

% Load coil reference scan
% If existing coil sensitivity matrices are saved in the GannetRecon_output
% folder, load them here. Otherwise, start loading them.
if (exist([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_sense_coil_img.mat'],'file') || ...
        exist([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_volume_coil_img.mat'],'file'))
    fprintf('Found existing coil reference data. Loading...\n%s\n',[spec_path filesep 'GannetRecon_output' filesep 'ref_scan_sense_coil_img.mat']);
    load([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_sense_coil_img.mat']);
    load([spec_path filesep 'GannetRecon_output' filesep 'ref_scan_volume_coil_img.mat']);
    disp('Loading coil reference data finished!');
else
    [ref_file,Ref,img_array,img_body,noise_array,noise_body] = loadRefScan(spec_path);
end

%% Calculate the actual unfolding matrix
MRS_struct = calcUnfoldingMatrix(MRS_struct, ref_file, Ref, img_array, noise_array, spec_path);

%% perform SENSE unfolding
% m=1 %1 for metab, 2 for water ref

% disp('sense unfolding...');
% 
% % clear signalunf signalres
% for m=1:2
%     for a =1:naveragesw1
%         signal = FID(:,:,m,a);
%         signalunf(:,:,m,a) = U*signal;
%         signalres(:,:,m,a) = signal - S*signalunf(:,:,m);
%     end
% end
    
% % Save all relevant data/information to MRS_struct
% MRS_struct.p.NVoxels=size(signalunf,1);
% MRS_struct.p.Navg = size(signalunf,4); % GO 11/01/2016
end
