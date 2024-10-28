function MRS_struct = GannetMask_GE_nii(fname, nii_file, MRS_struct, ii, vox, kk)

% Co-register GE P-files to structural images in NIfTI format. Code heavily
% based on Dr. Peter Van Schuerbbeek's (UZ Brussel) coreg_p code and Ralph
% Noeske's (GE Berlin) SV_MRI voxel co-registration code.

% Parse P-file to extract voxel geometry
if MRS_struct.p.GE.rdbm_rev_num(ii) >= 11.0
    fid = fopen(fname, 'r', 'ieee-le');
else
    fid = fopen(fname, 'r', 'ieee-be');
end

switch num2str(MRS_struct.p.GE.rdbm_rev_num(ii))
    case '14.3'
        rdb_hdr_off_image   = 377;
        rdb_hdr_ps_mps_freq = 107;
        tlhc                = 121;
        trhc                = 124;
        brhc                = 127;
    case '16'
        rdb_hdr_off_image   = 377;
        rdb_hdr_ps_mps_freq = 107;
        tlhc                = 133;
        trhc                = 136;
        brhc                = 139;
    case {'20.006','20.007','24'}
        rdb_hdr_off_image   = 377;
        rdb_hdr_ps_mps_freq = 107;
        tlhc                = 181;
        trhc                = 184;
        brhc                = 187;
    case {'26.002','27','27.001','28.002','28.003','30','30.1'}
        rdb_hdr_off_image   = 11;
        rdb_hdr_ps_mps_freq = 123;
        tlhc                = 181;
        trhc                = 184;
        brhc                = 187;
end

fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
o_hdr_value = fread(fid, brhc+2, 'real*4');
fclose(fid);

tlhc_RAS = o_hdr_value(tlhc:tlhc+2)';
trhc_RAS = o_hdr_value(trhc:trhc+2)';
brhc_RAS = o_hdr_value(brhc:brhc+2)';

e1_SVS_n = trhc_RAS - tlhc_RAS;
e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
e2_SVS_n = brhc_RAS - trhc_RAS;
e2_SVS_n = e2_SVS_n ./ norm(e2_SVS_n);
e3_SVS_n = -cross(e1_SVS_n, e2_SVS_n);

[~,orientation_SVS] = max(abs(e3_SVS_n));

if orientation_SVS == 3 % axial
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e2_SVS_n;
    e3_SVS_n2 = e3_SVS_n;
elseif orientation_SVS == 2 % coronal
    e1_SVS_n2 = e1_SVS_n;
    e2_SVS_n2 = e3_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
elseif orientation_SVS == 1 % sagittal
    e1_SVS_n2 = e3_SVS_n;
    e2_SVS_n2 = e1_SVS_n;
    e3_SVS_n2 = e2_SVS_n;
end

MRS_struct.p.voxang(ii,:) = get_euler(e1_SVS_n2, e2_SVS_n2, e3_SVS_n2);

e1_SVS = MRS_struct.p.voxdim(ii,1) * e1_SVS_n2;
e2_SVS = MRS_struct.p.voxdim(ii,2) * e2_SVS_n2;
e3_SVS = MRS_struct.p.voxdim(ii,3) * e3_SVS_n2;

% LPS gives center of voxel
LPS_SVS_edge = MRS_struct.p.voxoff(ii,:) - 0.5 * e1_SVS ...
                                         - 0.5 * e2_SVS ...
                                         - 0.5 * e3_SVS;

% Load in NIfTI file
V         = spm_vol(nii_file);
[T1, XYZ] = spm_read_vols(V);

% Shift imaging voxel coordinates by half an imaging voxel so that the XYZ matrix
% tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
[~,voxdim]      = spm_get_bbox(V,'fv');
voxdim          = abs(voxdim)';
halfpixshift    = -voxdim(1:3)/2;
halfpixshift(3) = -halfpixshift(3);
XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

dXYZ = sqrt((XYZ(1,:) - LPS_SVS_edge(1)).^2 + ...
            (XYZ(2,:) - LPS_SVS_edge(2)).^2 + ...
            (XYZ(3,:) - LPS_SVS_edge(3)).^2);
[~,refvox] = min(dXYZ);
[refvox_x, refvox_y, refvox_z] = ind2sub(V.dim, refvox(1));

e1_edge = LPS_SVS_edge + e1_SVS;
e2_edge = LPS_SVS_edge + e2_SVS;
e3_edge = LPS_SVS_edge + e3_SVS;

de1_XYZ = sqrt((XYZ(1,:) - e1_edge(1)).^2 + ...
               (XYZ(2,:) - e1_edge(2)).^2 + ...
               (XYZ(3,:) - e1_edge(3)).^2);
[~,e1vox] = min(de1_XYZ);
[e1vox_x, e1vox_y, e1vox_z] = ind2sub(V.dim, e1vox(1));

de2_XYZ = sqrt((XYZ(1,:) - e2_edge(1)).^2 + ...
               (XYZ(2,:) - e2_edge(2)).^2 + ...
               (XYZ(3,:) - e2_edge(3)).^2);
[~,e2vox] = min(de2_XYZ);
[e2vox_x, e2vox_y, e2vox_z] = ind2sub(V.dim, e2vox(1));

de3_XYZ = sqrt((XYZ(1,:) - e3_edge(1)).^2 + ...
               (XYZ(2,:) - e3_edge(2)).^2 + ...
               (XYZ(3,:) - e3_edge(3)).^2);
[~,e3vox] = min(de3_XYZ);
[e3vox_x, e3vox_y, e3vox_z] = ind2sub(V.dim, e3vox(1));

% Create a mask with all voxels that are inside the voxel
mask = zeros(V.dim);

nx = floor(sqrt((e1vox_x - refvox_x)^2 + (e1vox_y - refvox_y)^2 + (e1vox_z - refvox_z)^2)) * 2;
ny = floor(sqrt((e2vox_x - refvox_x)^2 + (e2vox_y - refvox_y)^2 + (e2vox_z - refvox_z)^2)) * 2;
nz = floor(sqrt((e3vox_x - refvox_x)^2 + (e3vox_y - refvox_y)^2 + (e3vox_z - refvox_z)^2)) * 2;

stepx = ([e1vox_x, e1vox_y, e1vox_z] - [refvox_x, refvox_y, refvox_z]) / nx;
stepy = ([e2vox_x, e2vox_y, e2vox_z] - [refvox_x, refvox_y, refvox_z]) / ny;
stepz = ([e3vox_x, e3vox_y, e3vox_z] - [refvox_x, refvox_y, refvox_z]) / nz;

mrs_box_ind = 1:(nx * ny * nz);
mrs_box_sub = zeros(3, nx * ny * nz);

[mrs_box_sub_x, mrs_box_sub_y, mrs_box_sub_z] = ind2sub([nx, ny, nz], mrs_box_ind);

mrs_box_sub(1,:) = mrs_box_sub_x;
mrs_box_sub(2,:) = mrs_box_sub_y;
mrs_box_sub(3,:) = mrs_box_sub_z;

e1_stepx = repmat(stepx, [numel(mrs_box_sub(1,:)), 1])';
e2_stepy = repmat(stepy, [numel(mrs_box_sub(1,:)), 1])';
e3_stepz = repmat(stepz, [numel(mrs_box_sub(1,:)), 1])';

mrs_box_sub = repmat((mrs_box_sub(1,:) - 1), [3 1]) .* e1_stepx + ...
              repmat((mrs_box_sub(2,:) - 1), [3 1]) .* e2_stepy + ...
              repmat((mrs_box_sub(3,:) - 1), [3 1]) .* e3_stepz;
refvox_rep  = repmat([refvox_x, refvox_y, refvox_z], [numel(mrs_box_sub(1,:)), 1])';
mrs_box_sub = round(mrs_box_sub + refvox_rep);
mrs_box_ind = sub2ind(V.dim, mrs_box_sub(1,:), ...
                             mrs_box_sub(2,:), ...
                             mrs_box_sub(3,:));
mask(mrs_box_ind) = 1;

% Build output (code to make voxel mask yellow borrowed from SPM12)

[a,b] = fileparts(fname);
if isempty(a)
    a = '.';
end
V_mask.fname   = fullfile([a filesep b '_mask.nii']);
V_mask.descrip = 'MRS_voxel_mask';
V_mask.dim     = V.dim;
V_mask.dt      = V.dt;
V_mask.mat     = V.mat;
V_mask         = spm_write_vol(V_mask, mask);

MRS_struct.mask.(vox{kk}).outfile(ii,:) = cellstr(V_mask.fname);

% Transform structural image and co-registered voxel mask from voxel to
% world space for output
voxel_ctr                = MRS_struct.p.voxoff(ii,:);
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


function euler_angles = get_euler(r1, r2, r3)

r1    = -r1;
r1(3) = -r1(3);
r2(3) = -r2(3);
r3(3) = -r3(3);

if abs(r3(1)) ~= 1
    theta1 = -asin(r3(1));
    %theta2 = pi - theta1;
    psi1 = atan2(r3(2)/cos(theta1), r3(3)/cos(theta1));
    %psi2 = atan2(r3(2)/cos(theta2), r3(3)/cos(theta2));
    phi1 = atan2(r2(1)/cos(theta1), r1(1)/cos(theta1));
    %phi2 = atan2(r2(1)/cos(theta2), r1(1)/cos(theta2));
else
    phi1 = 0;
    if r3(1) == -1
        theta1 = pi/2;
        psi1 = phi1 + atan2(r1(2), r1(3));
    else
        theta1 = -pi/2;
        psi1 = -phi1 + atan2(-r1(2), -r1(3));
    end
end

euler_angles(1) = round(-phi1*180/pi);
euler_angles(2) = round(-psi1*180/pi);
euler_angles(3) = round(theta1*180/pi);

end



