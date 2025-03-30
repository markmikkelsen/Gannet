function MRS_struct = NormalizeVoxelMask(MRS_struct, vox, ii, kk)
% Function to normalize binary MRS voxel masks to MNI space
%
% Author: Mark Mikkelsen, Ph.D. (Weill Cornell Medicine) (2024)

% Select forward deformation field(s) and MRS voxel mask(s) to be transformed into MNI space
voxel_mask = MRS_struct.mask.(vox{kk}).outfile{ii};
[vox_dir, vox_name, vox_ext] = fileparts(voxel_mask);
MRS_struct.mask.(vox{kk}).outfile_norm(ii,:) = cellstr(fullfile(vox_dir, ['w_' vox_name vox_ext]));

[struc_dir, struc_name, struc_ext] = fileparts(MRS_struct.mask.(vox{kk}).T1image{ii});
fwd_def = fullfile(struc_dir, ['y_' struc_name struc_ext]);

% Normalize MRS voxel mask to MNI space
matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr(fwd_def);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(voxel_mask);
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4; % 4th deg B-spline interpolation
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';

% Run job
spm_jobman('run',matlabbatch);



