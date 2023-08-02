function MRS_struct = GannetMask_GE(fname, dcm_dir, MRS_struct, ii, vox, kk)

% Updated GannetMask_GE function. Rotated localizer files are
% no longer required. Voxel geometry is taken directly from P-file headers
% and the structural image DICOMs.
% Code heavily based on Ralph Noeske's (GE Berlin) SV_MRI voxel
% co-registration code.

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
        image_user8         = 38;
        image_user11        = 41;
        tlhc                = 121;
        trhc                = 124;
        brhc                = 127;
    case '16'
        rdb_hdr_off_image   = 377;
        rdb_hdr_ps_mps_freq = 107;
        image_user8         = 50;
        image_user11        = 53;
        tlhc                = 133;
        trhc                = 136;
        brhc                = 139;
    case {'20.006','20.007','24'}
        rdb_hdr_off_image   = 377;
        rdb_hdr_ps_mps_freq = 107;
        image_user8         = 98;
        image_user11        = 101;
        tlhc                = 181;
        trhc                = 184;
        brhc                = 187;
    case {'26.002','27','27.001','28.002','28.003','30'}
        rdb_hdr_off_image   = 11;
        rdb_hdr_ps_mps_freq = 123;
        image_user8         = 98;
        image_user11        = 101;
        tlhc                = 181;
        trhc                = 184;
        brhc                = 187;
end

fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
o_hdr_value = fread(fid, brhc+2, 'real*4');
fclose(fid);

MRS_struct.p.voxdim(ii,:) = o_hdr_value(image_user8:image_user8+2)';
MRS_struct.p.voxoff(ii,:) = o_hdr_value(image_user11:image_user11+2)';
tlhc_RAS = o_hdr_value(tlhc:tlhc+2)';
trhc_RAS = o_hdr_value(trhc:trhc+2)';
brhc_RAS = o_hdr_value(brhc:brhc+2)';

% Convert from RAS to LPS
MRS_struct.p.voxoff(ii,:) = MRS_struct.p.voxoff(ii,:) .* [-1 -1 1];
tlhc_LPS = tlhc_RAS .* [-1 -1 1];
trhc_LPS = trhc_RAS .* [-1 -1 1];
brhc_LPS = brhc_RAS .* [-1 -1 1];

e1_SVS_n = trhc_LPS - tlhc_LPS;
e1_SVS_n = e1_SVS_n ./ norm(e1_SVS_n);
e2_SVS_n = brhc_LPS - trhc_LPS;
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

% Read all DICOM files into one volume
curr_dir = pwd;
data_dir = fileparts(MRS_struct.metabfile{1,ii});
if isempty(data_dir)
    data_dir = '.';
end
cd(data_dir);
data_dir = pwd;
cd(curr_dir);

if exist(dcm_dir, 'dir')
    cd(dcm_dir);
    dcm_list = dir;
    dcm_list = dcm_list(~ismember({dcm_list.name}, {'.','..','.DS_Store'}));
    dcm_list = cellstr(char(dcm_list.name));
    dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.nii'))); %#ok<*STRCLFH>
    dcm_list = dcm_list(cellfun(@isempty, strfind(dcm_list, '.mat')));
    dcm_hdr  = spm_dicom_headers(char(dcm_list));
    
    slice_location = zeros(1,length(dcm_list));
    for jj = 1:length(dcm_list)
        slice_location(jj) = dcm_hdr{jj}.SliceLocation;
    end
    
    % Order slices according to slice position
    [~,order_index] = sort(slice_location);
    tmp_hdr = dcm_hdr;
    for jj = 1:length(dcm_list)
        dcm_hdr{jj} = tmp_hdr{order_index(jj)};
    end
    
    cd(curr_dir);
elseif exist(dcm_dir, 'file')
    dcm_hdr = spm_dicom_headers(dcm_dir);
end

% Create NIfTI file of T1 image
nii_file_dir = spm_dicom_convert(dcm_hdr, 'all', 'flat', 'nii', data_dir);
nii_file     = nii_file_dir.files{1};
V            = spm_vol(nii_file);

if exist(dcm_dir, 'dir')
    MRI_voxel_size = [dcm_hdr{1}.PixelSpacing(1) ...
                      dcm_hdr{1}.PixelSpacing(2) ...
                      dcm_hdr{1}.SpacingBetweenSlices];
    
    MRI_dim = [dcm_hdr{1}.Rows ...
               dcm_hdr{1}.Columns ...
               dcm_hdr{1}.ImagesInAcquisition];
    
    e1_MRI_n = dcm_hdr{1}.ImageOrientationPatient(1:3);
    e2_MRI_n = dcm_hdr{1}.ImageOrientationPatient(4:6);
elseif exist(dcm_dir, 'file')
    MRI_voxel_size = [dcm_hdr{1}.SharedFunctionalGroupsSequence{1}.PixelMeasuresSequence{1}.PixelSpacing(1) ...
                      dcm_hdr{1}.SharedFunctionalGroupsSequence{1}.PixelMeasuresSequence{1}.PixelSpacing(2) ...
                      dcm_hdr{1}.SharedFunctionalGroupsSequence{1}.PixelMeasuresSequence{1}.SliceThickness];
    
    MRI_dim = [dcm_hdr{1}.Rows ...
               dcm_hdr{1}.Columns ...
               dcm_hdr{1}.NumberOfFrames];
    
    e1_MRI_n = dcm_hdr{1}.PerFrameFunctionalGroupsSequence{1}.PlaneOrientationSequence{1}.ImageOrientationPatient(1:3);
    e2_MRI_n = dcm_hdr{1}.PerFrameFunctionalGroupsSequence{1}.PlaneOrientationSequence{1}.ImageOrientationPatient(4:6);
end
e3_MRI_n = cross(e1_MRI_n, e2_MRI_n); % e3 vector is perpendicular to the slice orientation

[~,orientation_MRI] = max(abs(e3_MRI_n));
if orientation_MRI == 2 % coronal
    e3_MRI_n = -e3_MRI_n;
end

if orientation_MRI == 3 % axial
    e1_MRI_n2      = e1_MRI_n;
    e2_MRI_n2      = e2_MRI_n;
    e3_MRI_n2      = e3_MRI_n;
    MRI_voxel_size = MRI_voxel_size([2 1 3]);
    MRI_dim        = MRI_dim([2 1 3]);
elseif orientation_MRI == 2 % coronal
    e1_MRI_n2      = e1_MRI_n;
    e2_MRI_n2      = e3_MRI_n;
    e3_MRI_n2      = e2_MRI_n;
    MRI_voxel_size = MRI_voxel_size([2 3 1]);
    MRI_dim        = MRI_dim([2 3 1]);
elseif orientation_MRI == 1 % sagittal
    e1_MRI_n2      = e3_MRI_n;
    e2_MRI_n2      = e1_MRI_n;
    e3_MRI_n2      = e2_MRI_n;
    MRI_voxel_size = MRI_voxel_size([3 2 1]);
    MRI_dim        = MRI_dim([3 2 1]);
end

% LPS_edge gives location of the edge of the image volume
if exist(dcm_dir, 'dir')
    LPS_MRI_center = dcm_hdr{1}.ImagePositionPatient;
elseif exist(dcm_dir, 'file')
    LPS_MRI_center = dcm_hdr{1}.PerFrameFunctionalGroupsSequence{1}.PlanePositionSequence{1}.ImagePositionPatient;
end
LPS_MRI_edge   = LPS_MRI_center - 0.5 * MRI_voxel_size(1) * e1_MRI_n2 ...
                                - 0.5 * MRI_voxel_size(2) * e2_MRI_n2 ...
                                - 0.5 * MRI_voxel_size(3) * e3_MRI_n2;

% Create voxel mask
E_MRI = [e1_MRI_n2 e2_MRI_n2 e3_MRI_n2];
c_MRS = MRS_struct.p.voxoff(ii,:)';
c_MRI = E_MRI' * (c_MRS - LPS_MRI_edge);
d_MRI = c_MRI ./ MRI_voxel_size';
s_MRS = sqrt(sum(MRS_struct.p.voxdim(ii,:).^2))/2;
d_MRS = s_MRS ./ MRI_voxel_size';

[Xm,Ym,Zm] = ndgrid(1:MRI_dim(1), 1:MRI_dim(2), 1:MRI_dim(3));
X = LPS_MRI_center(1) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(1) ...
                      + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(1) ...
                      + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(1);
Y = LPS_MRI_center(2) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(2) ...
                      + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(2) ...
                      + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(2);
Z = LPS_MRI_center(3) + (Xm-1) * MRI_voxel_size(1) * e1_MRI_n2(3) ...
                      + (Ym-1) * MRI_voxel_size(2) * e2_MRI_n2(3) ...
                      + (Zm-1) * MRI_voxel_size(3) * e3_MRI_n2(3);

P_1  = LPS_SVS_edge;
P_2  = LPS_SVS_edge + e1_SVS; % L
P_3  = LPS_SVS_edge + e2_SVS; % P
P_4  = LPS_SVS_edge + e3_SVS; % S
A    = zeros(3,1);
mask = zeros(MRI_dim);

for e1 = max(floor(d_MRI(1) - d_MRS(1)), 0):min(ceil(d_MRI(1) + d_MRS(1)), size(mask,1)) % L
    for e2 = max(floor(d_MRI(2) - d_MRS(2)), 0):min(ceil(d_MRI(2) + d_MRS(2)), size(mask,2)) % P
        for e3 = max(floor(d_MRI(3) - d_MRS(3)), 0):min(ceil(d_MRI(3) + d_MRS(3)), size(mask,3)) % S
            A(1) = X(e1,e2,e3);
            A(2) = Y(e1,e2,e3);
            A(3) = Z(e1,e2,e3);
            % Distance of A to planes in SI direction
            d_5 = e3_SVS_n2 * (A - P_1');
            d_6 = -e3_SVS_n2 * (A - P_4');
            if d_5 >= 0 && d_6 >= 0
                % Distance of A to planes in AP direction
                d_3 = e2_SVS_n2 * (A - P_1');
                d_4 = -e2_SVS_n2 * (A - P_3');
                if d_3 >= 0 && d_4 >= 0
                    % Distance of A to planes in RL direction
                    d_1 = e1_SVS_n2 * (A - P_1');
                    d_2 = -e1_SVS_n2 * (A - P_2');
                    if d_1 >= 0 && d_2 >= 0
                        mask(e1,e2,e3) = 1;
                    end
                end
            end
        end
    end
end

% Convert from LPS to RAS
if orientation_MRI == 2 % coronal
    mask = permute(mask, [1 3 2]);
    mask = flip(mask,3);
elseif orientation_MRI == 1 % sagittal
    mask = permute(mask, [2 3 1]);
end
mask = flip(mask,2);

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
voxel_ctr(1:2)           = -voxel_ctr(1:2);
[img_t, img_c, img_s]    = voxel2world_space(V, voxel_ctr);
[mask_t, mask_c, mask_s] = voxel2world_space(V_mask, voxel_ctr);

T1 = spm_read_vols(V);

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



