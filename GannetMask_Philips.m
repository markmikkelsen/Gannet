function MRS_struct = GannetMask_Philips(fname, nii_file, MRS_struct, ii, vox, kk)

if nargin == 2
    MRS_struct.ii = 1;
    ii = 1;
end

[path, name] = fileparts(fname);
fidoutmask   = fullfile(path, [name '_mask.nii']);

[~,~,ext] = fileparts(fname);
if all(isstrprop(ext(end-3:end), 'upper'))
    spar_ext = 'SPAR';
else
    spar_ext = 'spar';
end
sparname     = [fname(1:(end-4)) spar_ext];
sparname     = fopen(sparname,'r');
sparheadinfo = textscan(sparname, '%s');
sparheadinfo = sparheadinfo{1};

sparidx = find(ismember(sparheadinfo, 'ap_size') == 1);
MRS_struct.p.voxdim(ii,2) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'lr_size') == 1);
MRS_struct.p.voxdim(ii,1) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'cc_size') == 1);
MRS_struct.p.voxdim(ii,3) = str2double(sparheadinfo{sparidx + 2});

sparidx = find(ismember(sparheadinfo, 'ap_off_center') == 1);
MRS_struct.p.voxoff(ii,2) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'lr_off_center') == 1);
MRS_struct.p.voxoff(ii,1) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'cc_off_center') == 1);
MRS_struct.p.voxoff(ii,3) = str2double(sparheadinfo{sparidx + 2});

sparidx = find(ismember(sparheadinfo, 'ap_angulation') == 1);
MRS_struct.p.voxang(ii,2) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'lr_angulation') == 1);
MRS_struct.p.voxang(ii,1) = str2double(sparheadinfo{sparidx + 2});
sparidx = find(ismember(sparheadinfo, 'cc_angulation') == 1);
MRS_struct.p.voxang(ii,3) = str2double(sparheadinfo{sparidx + 2});

V         = spm_vol(nii_file);
[T1, XYZ] = spm_read_vols(V);

% Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
% tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim]      = spm_get_bbox(V,'fv');
voxdim          = abs(voxdim)';
halfpixshift    = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ             = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

% Get information from SPAR - change later to be read in
ap_size = MRS_struct.p.voxdim(ii,2);
lr_size = MRS_struct.p.voxdim(ii,1);
cc_size = MRS_struct.p.voxdim(ii,3);
ap_off  = MRS_struct.p.voxoff(ii,2);
lr_off  = MRS_struct.p.voxoff(ii,1);
cc_off  = MRS_struct.p.voxoff(ii,3);
ap_ang  = MRS_struct.p.voxang(ii,2);
lr_ang  = MRS_struct.p.voxang(ii,1);
cc_ang  = MRS_struct.p.voxang(ii,3);

% We need to flip ap and lr axes to match NIfTI convention
ap_off = -ap_off;
lr_off = -lr_off;
ap_ang = -ap_ang;
lr_ang = -lr_ang;

% Define the voxel - use x y z
% Currently have spar convention that have in AUD voxel - will need to
% check for everything in future...
% x - left = positive
% y - posterior = postive
% z - superior = positive
vox_ctr = ...
    [lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2 -ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2  cc_size/2;
     lr_size/2  ap_size/2  cc_size/2;
    -lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2  ap_size/2 -cc_size/2;
     lr_size/2 -ap_size/2 -cc_size/2;
    -lr_size/2 -ap_size/2 -cc_size/2];

% Make rotations on voxel
rad = pi/180;

xrot      = zeros(3,3);
xrot(1,1) = 1;
xrot(2,2) = cos(lr_ang * rad);
xrot(2,3) = -sin(lr_ang * rad);
xrot(3,2) = sin(lr_ang * rad);
xrot(3,3) = cos(lr_ang * rad);

yrot      = zeros(3,3);
yrot(1,1) = cos(ap_ang * rad);
yrot(1,3) = sin(ap_ang * rad);
yrot(2,2) = 1;
yrot(3,1) = -sin(ap_ang * rad);
yrot(3,3) = cos(ap_ang * rad);

zrot      = zeros(3,3);
zrot(1,1) = cos(cc_ang * rad);
zrot(1,2) = -sin(cc_ang * rad);
zrot(2,1) = sin(cc_ang * rad);
zrot(2,2) = cos(cc_ang * rad);
zrot(3,3) = 1;

% Rotate voxel
vox_rot = xrot * yrot * zrot * vox_ctr.';

% Calculate corner coordinates relative to xyz origin
vox_ctr_coor = [lr_off ap_off cc_off];
vox_ctr_coor = repmat(vox_ctr_coor.', [1,8]);
vox_corner   = vox_rot + vox_ctr_coor;

mask          = zeros(1,size(XYZ,2));
sphere_radius = sqrt((lr_size/2).^2 + (ap_size/2).^2 + (cc_size/2).^2);
dist2voxctr   = sqrt(sum((XYZ - repmat([lr_off ap_off cc_off].', [1 size(XYZ,2)])).^2, 1));
sphere_mask(dist2voxctr <= sphere_radius) = 1;

mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);

tri      = delaunayn([vox_corner.'; [lr_off ap_off cc_off]]);
tn       = tsearchn([vox_corner.'; [lr_off ap_off cc_off]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask == 1) = isinside;

mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;
V_mask         = spm_write_vol(V_mask, mask);

% Build output (code to make voxel mask yellow borrowed from SPM12)

fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.(vox{kk}).outfile(ii,:) = fidoutmask;

% Transform structural image and co-registered voxel mask from voxel to
% world space for output
voxel_ctr                = [lr_off ap_off cc_off];
[img_t, img_c, img_s]    = voxel2world_space(V, voxel_ctr);
[mask_t, mask_c, mask_s] = voxel2world_space(V_mask, voxel_ctr);

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



