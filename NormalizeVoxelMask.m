function MRS_struct = NormalizeVoxelMask(MRS_struct, vox, ii, kk)
% Function to normalize binary MRS voxel masks to MNI space
%
% Author: Mark Mikkelsen, Ph.D. (Weill Cornell Medicine) (2024)

mergestructs = @(x,y) cell2struct([struct2cell(x); struct2cell(y)], [fieldnames(x); fieldnames(y)]);

% Select forward deformation field(s) and MRS voxel mask(s) to be transformed into MNI space
vox_mask = MRS_struct.mask.(vox{kk}).fname{ii};
[vox_dir, vox_name, vox_ext] = fileparts(vox_mask);
MRS_struct.mask.(vox{kk}).fname_norm(ii,:) = cellstr(fullfile(vox_dir, ['w_' vox_name vox_ext]));

fwd_def = MRS_struct.mask.(vox{kk}).fwd_def{ii};

% Normalize MRS voxel mask to MNI space
matlabbatch{1}.spm.spatial.normalise.write.subj.def        = cellstr(fwd_def);
matlabbatch{1}.spm.spatial.normalise.write.subj.resample   = cellstr(vox_mask);
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox    = [2 2 2];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4; % 4th deg B-spline interpolation
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w_';

% Run job
spm_jobman('run',matlabbatch);

% BIDSify
if MRS_struct.p.bids
    bids_file = bids.File(vox_mask);
    input = mergestructs(bids_file.entities, struct('space', 'MNI152'));
    bids_file.entities = input;
    vox_mask_norm = fullfile(MRS_struct.out.BIDS.pth, 'derivatives', 'Gannet_output', bids_file.bids_path, bids_file.filename);
    movefile(MRS_struct.mask.(vox{kk}).fname_norm{ii}, vox_mask_norm);
    MRS_struct.mask.(vox{kk}).fname_norm{ii} = vox_mask_norm;
end



