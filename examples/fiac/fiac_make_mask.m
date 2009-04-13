function V_rmsk = fiac_make_mask(fname, defaults)
% Make mask for within-brain voxels from SPM gray/white/csf images
% and reslice into space of default normalized output image

sdir = fullfile(spm('dir'), 'apriori');
P1 = spm_get('files', sdir, 'gray.mnc');
P2 = spm_get('files', sdir, 'white.mnc');
P3 = spm_get('files', sdir, 'csf.mnc');
V = spm_vol(strvcat(P1, P2, P3));

img = zeros(V(1).dim(1:3));
for i_no = 1:length(V)
  img = img + spm_read_vols(V(i_no));
end
V_msk = V(1);
V_msk.fname = fname;
V_msk = spm_write_vol(V_msk, img);

nd = defaults.normalise.write;
[out_dim out_mat] = sf_bb_to_dim_mat(V_msk.mat, nd.bb, nd.vox);
M = V_msk.mat \ out_mat;
img = zeros(out_dim);
for p = 1:out_dim(3)
  img(:,:,p) = spm_slice_vol(V_msk, M * spm_matrix([0 0 p]), out_dim(1:2), 0);
end
V_rmsk = V_msk;
V_rmsk.dim(1:3) = out_dim;
V_rmsk.mat = out_mat;
V_rmsk = spm_write_vol(V_rmsk, img);

return

function [dim, mat] = sf_bb_to_dim_mat(template_mat, bb,vox)
% Calculate dimensions, mat from template mat, bounding box, voxels
% FORMAT [dim, mat] = sf_bb_to_dim_mat(template_mat, bb,vox)
% 
% Stolen from get_xyzmat subfunction of spm_write_sn

msk       = find(vox<0);
bb        = sort(bb);
bb(:,msk) = flipud(bb(:,msk));

% Adjust bounding box slightly - so it rounds to closest voxel.
bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

M   = template_mat;
vxg = sqrt(sum(M(1:3,1:3).^2));
if det(M(1:3,1:3))<0, vxg(1) = -vxg(1); end;
ogn = M\[0 0 0 1]';
ogn = ogn(1:3)';

% Calculate dimensions
dim = floor(diff(bb) ./ vox) + 1;

% mat transform
og  = -vxg.*ogn;
of  = -vox.*(round(-bb(1,:)./vox)+1);
M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
mat = template_mat*inv(M1)*M2;

if (spm_flip_analyze_images & det(mat(1:3,1:3))>0) | (~spm_flip_analyze_images & det(mat(1:3,1:3))<0),
	Flp = [-1 0 0 (length(x)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
	mat = mat*Flp;
end;
return;














