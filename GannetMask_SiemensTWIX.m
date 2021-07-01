function MRS_struct = GannetMask_SiemensTWIX(fname, nii_file, MRS_struct, ii, vox, kk)
%   Creates a .nii file containing the voxel mask of the MRS voxel.
%   Needs to be called from GannetCoRegister.
%   Requires SPM8 or SPM12 to be added to the MATLAB path.
%
%   Author:
%       Dr. Georg Oeltzschner (Johns Hopkins University, 2018-02-24)
%       goeltzs1@jhmi.edu
%   
%   Credits:
%       The routine for correct determination of the phase and readout
%       directions of the MRS voxel is adapted from 
%       vox2ras_rsolveAA.m
%       (Dr. Rudolph Pienaar, Massachusetts General Hospital, Boston)
%
%   History:
%       2018-02-24: New version of GannetMask_SiemensTWIX
%       2020-12-04: Updates to voxel dipaly in output figure

[path, name] = fileparts(fname);
fidoutmask   = fullfile(path,[name '_mask.nii']);

% Extract voxel position and rotation parameters from MRS_struct
NormSag = MRS_struct.p.NormSag(ii);
NormCor = MRS_struct.p.NormCor(ii);
NormTra = MRS_struct.p.NormTra(ii);
VoI_InPlaneRot = MRS_struct.p.VoI_InPlaneRot(ii);
% Correct voxel offsets by table position (if field exists)
if isfield(MRS_struct.p,'TablePosition')
    VoxOffs = [MRS_struct.p.voxoff(ii,1) + MRS_struct.p.TablePosition(ii,1) ...
               MRS_struct.p.voxoff(ii,2) + MRS_struct.p.TablePosition(ii,2) ...
               MRS_struct.p.voxoff(ii,3) + MRS_struct.p.TablePosition(ii,3)];
else
    VoxOffs = [MRS_struct.p.voxoff(ii,1) MRS_struct.p.voxoff(ii,2) MRS_struct.p.voxoff(ii,3)];
end


% Parse direction cosines of the MRS voxel's normal vector and the rotation angle
% around the normal vector
% The direction cosine is the cosine of the angle between the normal
% vector and the respective direction.
% Example: If the normal vector points exactly along the FH direction, then:
% NormSag = cos(90) = 0, NormCor = cos(90) = 0, NormTra = cos(0) = 1.
Norm = [-NormSag -NormCor NormTra];
ROT = VoI_InPlaneRot;
% Find largest element of normal vector of the voxel to determine primary
% orientation.
% Example: if NormTra has the smallest out of the three Norm
% values, the angle of the normal vector with the Tra direction (FH) is the
% smallest, and the primary orientation is transversal.
[~, maxdir] = max([abs(NormSag) abs(NormCor) abs(NormTra)]);
switch maxdir
    case 1
        vox_orient = 's'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 2
        vox_orient = 'c'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
    case 3
        vox_orient = 't'; % 't' = transversal, 's' = sagittal', 'c' = coronal;
end

% Phase reference vector
% Adapted from Rudolph Pienaar's "vox2ras_rsolveAA.m" and
% Andre van der Kouwe's "autoaligncorrect.cpp"
Phase = zeros(3,1);
switch vox_orient
    case 't'
        % For transversal voxel orientation, the phase reference vector lies in
        % the sagittal plane
        Phase(1) = 0;
        Phase(2) = Norm(3)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        Phase(3) = -Norm(2)*sqrt(1/(Norm(2)*Norm(2)+Norm(3)*Norm(3)));
        VoxDims  = [MRS_struct.p.voxdim(ii,1) MRS_struct.p.voxdim(ii,2) MRS_struct.p.voxdim(ii,3)];
    case 'c'
        % For coronal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1) = Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2) = -Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3) = 0;
        VoxDims  = [MRS_struct.p.voxdim(ii,1) MRS_struct.p.voxdim(ii,2) MRS_struct.p.voxdim(ii,3)];
    case 's'
        % For sagittal voxel orientation, the phase reference vector lies in
        % the transversal plane
        Phase(1) = -Norm(2)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(2) = Norm(1)*sqrt(1/(Norm(1)*Norm(1)+Norm(2)*Norm(2)));
        Phase(3) = 0;
        VoxDims  = [MRS_struct.p.voxdim(ii,1) MRS_struct.p.voxdim(ii,2) MRS_struct.p.voxdim(ii,3)];
end

% The readout reference vector is the cross product of Norm and Phase
Readout    = cross(Norm, Phase);
M_R        = zeros(4,4);
M_R(1:3,1) = Phase;
M_R(1:3,2) = Readout;
M_R(1:3,3) = Norm;

% Define matrix for rotation around in-plane rotation angle
M3_Mu = [cos(ROT) sin(ROT) 0
        -sin(ROT) cos(ROT) 0
         0        0        1];

M3_R = M_R(1:3,1:3) * M3_Mu;
M_R(1:3,1:3) = M3_R;

% The MGH vox2ras matrix inverts the Readout column
M_R	= M_R * [1  0  0  0
             0 -1  0  0
             0  0  1  0
             0  0  0  1];

% Final rotation matrix
rotmat = M_R(1:3,1:3);

V        = spm_vol(nii_file);
[T1,XYZ] = spm_read_vols(V);
MRS_struct.mask.(vox{kk}).T1max(ii) = max(T1(:));

% Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
% tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim] = spm_get_bbox(V,'fv');
voxdim = abs(voxdim)';
halfpixshift = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

% We need to flip ap and lr axes to match NIFTI convention
VoxOffs(1) = -VoxOffs(1);
VoxOffs(2) = -VoxOffs(2);

% Define voxel coordinates before rotation and transition
vox_ctr = ...
    [VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2  VoxDims(3)/2;
    -VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2  VoxDims(2)/2 -VoxDims(3)/2;
     VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2;
    -VoxDims(1)/2 -VoxDims(2)/2 -VoxDims(3)/2];

% Apply rotation as prescribed
vox_rot = rotmat * vox_ctr.';

% Shift rotated voxel by the center offset to its final position
vox_ctr_coor = [VoxOffs(1) VoxOffs(2) VoxOffs(3)];
vox_ctr_coor = repmat(vox_ctr_coor.', [1 8]);
vox_corner = vox_rot + vox_ctr_coor;

% Create a mask with all voxels that are inside the voxel
mask = zeros(1,size(XYZ,2));
sphere_radius = sqrt((VoxDims(1)/2)^2 + (VoxDims(2)/2)^2 + (VoxDims(3)/2)^2);
distance2voxctr = sqrt(sum((XYZ-repmat([VoxOffs(1) VoxOffs(2) VoxOffs(3)].',[1 size(XYZ, 2)])).^2,1));
sphere_mask(distance2voxctr <= sphere_radius) = 1;
mask(sphere_mask == 1) = 1;
XYZ_sphere = XYZ(:,sphere_mask == 1);
tri = delaunayn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]]);
tn = tsearchn([vox_corner.'; [VoxOffs(1) VoxOffs(2) VoxOffs(3)]], tri, XYZ_sphere.');
isinside = ~isnan(tn);
mask(sphere_mask == 1) = isinside;

% Take over the voxel dimensions from the structural
mask = reshape(mask, V.dim);

V_mask.fname   = fidoutmask ;
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;
V_mask         = spm_write_vol(V_mask, mask);

% Build output (code to make voxel mask yellow borrowed from SPM12)

fidoutmask = cellstr(fidoutmask);
MRS_struct.mask.(vox{kk}).outfile(ii,:) = fidoutmask;
% Not clear how to formulate the rotations for triple rotations (revisit)
MRS_struct.p.voxang(ii,:) = [NaN NaN NaN];

% Transform structural image and co-registered voxel mask from voxel to
% world space for output
[img_t, img_c, img_s]    = voxel2world_space(V, VoxOffs);
[mask_t, mask_c, mask_s] = voxel2world_space(V_mask, VoxOffs);

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



