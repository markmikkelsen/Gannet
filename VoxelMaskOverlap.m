function MRS_struct = VoxelMaskOverlap(MRS_struct)
% Function to create a voxel overlap mask in MNI space
%
% Code borrowed from Osprey function OspreySeg.m (https://github.com/schorschinho/osprey)

% CREDIT:
% Helge Zollner, Johns Hopkins University (2022)

if MRS_struct.p.PRIAM
    vox = MRS_struct.p.vox;
else
    vox = MRS_struct.p.vox(1);
end

for kk = 1:length(vox)

    if MRS_struct.p.append
        out_dir = fullfile(pwd, 'Gannet_output');
    else
        if ~exist(fullfile(pwd, 'GannetMask_output'), 'dir')
            mkdir(fullfile(pwd, 'GannetMask_output'));
        end
        out_dir = fullfile(pwd, 'GannetMask_output');
    end

    MRS_struct.mask.(vox{kk}).mask_overlap.filename = fullfile(out_dir, 'MRS_voxel_overlap.nii');

    % The expression to calculate the average over all subjects' binary voxel masks
    expression = '(';
    for jj = 1:MRS_struct.p.numScans
        expression = [expression 'i' num2str(jj)]; %#ok<*AGROW>
        if jj ~= MRS_struct.p.numScans
            expression = [expression '+'];
        else
            expression = [expression ')/' num2str(MRS_struct.p.numScans)];
        end
    end
    MRS_struct.mask.(vox{kk}).mask_overlap.expression = expression;

    % Calculate an average overlap mask from subjects' voxel masks
    matlabbatch{1}.spm.util.imcalc.input          = cellstr(MRS_struct.mask.(vox{kk}).outfile_norm);
    matlabbatch{1}.spm.util.imcalc.output         = 'MRS_voxel_overlap';
    matlabbatch{1}.spm.util.imcalc.outdir         = cellstr(out_dir);
    matlabbatch{1}.spm.util.imcalc.expression     = expression;
    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = -4; % 4th deg sinc interpolation
    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;

    % Run job
    spm_jobman('run',matlabbatch);

end



