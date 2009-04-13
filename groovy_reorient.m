function groovy_reorient(glob_ps, sub_ps)
% add some transformation for all subject images
%
% Transformations stored in sub_ps.trans
% 4x4 matrix -> set the orientation to this value
% translation, rotation etc vector -> add to current

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  trans = this_sub.trans;
  
  % Do nothing if trans == 0
  if prod(size(trans)) == 1
    if trans == 0, continue, end
  end
  
  % 4x4 matrix -> set the orientation to this value
  % translation, rotation etc vector -> add to current
  set_to_f = all(size(trans)==4);
  if ~set_to_f, trans = spm_matrix(trans); end
  
  r_filter = [glob_ps.prefixes.realign this_sub.raw_filter];
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_get('files', dirn, r_filter);
    for i = 1:size(P,1)
      fname = deblank(P(i,:));
      if set_to_f
	spm_get_space(fname, trans);
      else
	M = spm_get_space(fname);
	spm_get_space(fname, trans * M);
      end
    end
  end
  
end











