function MRS_struct = GannetMask_NIfTI(fname, nii_file, MRS_struct, ii, vox, kk)

% Co-register NIfTI-MRS files to structural images in NIfTI format. Code
% heavily based on coreg_nifti.m from Osprey.

% CREDITS:
% Chris Davies-Jenkins, Johns Hopkins University 2022.
% Xiangrui Li, Ph.D. for his helpful suggestions using nii_tool.

nii_struc  = nii_tool('load', nii_file); % load structural nifti
nii_mrsvox = nii_tool('load', fname);    % load voxel nifti

% nii_viewer(NiiStruct, NiiVox); % overlay voxel on structural

% Assume MRS voxel and structural are in same space
nii_mrsvox.hdr.sform_code = nii_struc.hdr.sform_code;
nii_mrsvox.hdr.qform_code = nii_struc.hdr.qform_code;

nii_mrsvox.img            = 1; % overwrites image, so mask
nii_mrsvox.hdr.dim(4:end) = 1; % remove additional MRS dimensions from header

[a,b,c] = fileparts(fname);
if isempty(a)
    a = '.';
end
if strcmpi(c,'.gz')
    b(end-3:end) = [];
end
mask_filename = fullfile([a filesep b '_mask.nii']);
% Transform voxel to image resolution and save under mask_filename
nii_xform(nii_mrsvox, nii_struc.hdr, mask_filename, 'linear', 0);

% Load structural using SPM
V  = spm_vol(nii_file);
T1 = spm_read_vols(V);
MRS_struct.mask.(vox{kk}).T1max(ii) = max(T1(:));

% Load mask using SPM to adapt some fields
V_mask         = spm_vol(mask_filename);
V_mask.dt      = V.dt;
V_mask.descrip = 'MRS_voxel_mask';
V_mask         = spm_write_vol(V_mask, V_mask.private.dat(:,:,:)); % write mask

MRS_struct.mask.(vox{kk}).outfile(ii,:) = cellstr(V_mask.fname);
% Not clear how to formulate the rotations for triple rotations (revisit)
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN];

voxel_ctr = [nii_mrsvox.hdr.qoffset_x, nii_mrsvox.hdr.qoffset_y, nii_mrsvox.hdr.qoffset_z]; % CWDJ need to verify this

MRS_struct.p.voxoff(ii,:) = voxel_ctr;
[img_t, img_c, img_s]     = voxel2world_space(V, voxel_ctr);
[mask_t, mask_c, mask_s]  = voxel2world_space(V_mask, voxel_ctr);

w_t = zeros(size(img_t));
w_c = zeros(size(img_c));
w_s = zeros(size(img_s));

T1    = T1(:);
img_t = repmat(img_t / (mean(T1(T1 > 0.01)) + 3*std(T1(T1 > 0.01))), [1 1 3]);
img_c = repmat(img_c / (mean(T1(T1 > 0.01)) + 3*std(T1(T1 > 0.01))), [1 1 3]);
img_s = repmat(img_s / (mean(T1(T1 > 0.01)) + 3*std(T1(T1 > 0.01))), [1 1 3]);

c_img_t = zeros(size(img_t));
c_img_c = zeros(size(img_c));
c_img_s = zeros(size(img_s));

vox_mx = 1;
vox_mn = 0;

mask_t(mask_t(:) < vox_mn) = vox_mn;
mask_t(mask_t(:) > vox_mx) = vox_mx;
mask_t = (mask_t - vox_mn) / (vox_mx - vox_mn);

mask_c(mask_c(:) < vox_mn) = vox_mn;
mask_c(mask_c(:) > vox_mx) = vox_mx;
mask_c = (mask_c - vox_mn) / (vox_mx - vox_mn);

mask_s(mask_s(:) < vox_mn) = vox_mn;
mask_s(mask_s(:) > vox_mx) = vox_mx;
mask_s = (mask_s - vox_mn) / (vox_mx - vox_mn);

mask_t = 0.4 * mask_t;
mask_c = 0.4 * mask_c;
mask_s = 0.4 * mask_s;

vox_color = [1 1 0];

c_img_t = c_img_t + cat(3, mask_t * vox_color(1,1), mask_t * vox_color(1,2), mask_t * vox_color(1,3));
c_img_c = c_img_c + cat(3, mask_c * vox_color(1,1), mask_c * vox_color(1,2), mask_c * vox_color(1,3));
c_img_s = c_img_s + cat(3, mask_s * vox_color(1,1), mask_s * vox_color(1,2), mask_s * vox_color(1,3));

w_t = w_t + mask_t;
w_c = w_c + mask_c;
w_s = w_s + mask_s;

img_t = repmat(1 - w_t, [1 1 3]) .* img_t + c_img_t;
img_c = repmat(1 - w_c, [1 1 3]) .* img_c + c_img_c;
img_s = repmat(1 - w_s, [1 1 3]) .* img_s + c_img_s;

img_t(img_t < 0) = 0; img_t(img_t > 1) = 1;
img_c(img_c < 0) = 0; img_c(img_c > 1) = 1;
img_s(img_s < 0) = 0; img_s(img_s > 1) = 1;

img_t = flipud(img_t);
img_c = flipud(img_c);
img_s = flipud(img_s);

size_max        = max([max(size(img_t)) max(size(img_c)) max(size(img_s))]);
three_plane_img = zeros([size_max 3*size_max 3]);
three_plane_img(:,1:size_max,:)              = image_center(img_t, size_max);
three_plane_img(:,size_max+(1:size_max),:)   = image_center(img_s, size_max);
three_plane_img(:,size_max*2+(1:size_max),:) = image_center(img_c, size_max);

MRS_struct.mask.(vox{kk}).img{ii}       = three_plane_img;
MRS_struct.mask.(vox{kk}).T1image(ii,:) = {nii_file};

end



