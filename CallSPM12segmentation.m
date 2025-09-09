function CallSPM12segmentation(anatimage)

spmhome = fileparts(which('spm'));
tpm     = cellstr(spm_select('ExtFPList',fullfile(spmhome,'tpm'),'TPM.nii'));

matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[anatimage ',1']};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm    = tpm(1);
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm    = tpm(2);
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm    = tpm(3);
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm    = tpm(4);
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus  = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm    = tpm(5);
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus  = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm    = tpm(6);
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus  = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write       = [0 1];

spm_jobman('run',matlabbatch);

end



