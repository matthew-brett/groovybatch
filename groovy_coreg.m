function groovy_coreg(glob_ps, sub_ps)
% metabatch file for coregistration
  
flags = glob_ps.defaults.coreg;

for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  PG = this_sub.coreg.target;
  VG = spm_vol(PG);

  PF = this_sub.coreg.object;
  VF = spm_vol(PF);
   
  % other images (e.g., coregister the grey matter segment)
  PO = PF;
  for od = 1:length(this_sub.coreg.other_dirs)
    PO = strvcat(PO, ...
		 spm_get('files', this_sub.coreg.other_dirs{od}, ...
			 this_sub.coreg.other_filter));
  end
  PO = char(unique(cellstr(PO)));
  
  % do coregistration
  if ~isfield(glob_ps, 'coreg')
    glob_ps.coreg.method='spm';
  end
  if ~isfield(glob_ps.coreg, 'method')
    glob_ps.coreg.method='spm';
  end
  switch glob_ps.coreg.method
   case 'spm'
    x  = spm_coreg(VG, VF, flags.estimate);
    M  = inv(spm_matrix(x));
   case 'flirt'
    A = flirt(VF, VG, 'normcorr');
    M = inv(VF.mat*A/VG.mat);
   case 'flirt_9'
    A = flirt(VF, VG, 'normcorr', '-dof 9');
    M = inv(VF.mat*A/VG.mat);
   otherwise
    error('Coregistration is a bit skew-whiff');
  end
  MM = zeros(4,4,size(PO,1));
  for j=1:size(PO,1),
    MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
  end
  for j=1:size(PO,1),
    spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
  end
end
