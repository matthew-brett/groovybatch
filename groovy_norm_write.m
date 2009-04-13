function groovy_norm_write(glob_ps, sub_ps)
% metabatch file to write with normalization parameters for SPM2

defs = glob_ps.defaults.normalise;

for s = 1:length(sub_ps) % for each subject 
  my_sub = sub_ps(s);
  subj_dir = fullfile(glob_ps.fdata_root, my_sub.dir);
  
  % Make the default normalization parameters file name
  ns = my_sub.norm.source; 
  matname = [spm_str_manip(ns, 'sd') '_sn.mat'];
  
  % Do the reslicing for miscellaneous images
  spm_write_sn(my_sub.norm.others,matname,defs.write);

  % Do the reslicing for the EPIs
  nw_filter = [glob_ps.prefixes.norm_write my_sub.raw_filter];
  imgs = '';
  for ss = 1:length(my_sub.sesses)
    dirn = fullfile(subj_dir,my_sub.sesses(ss).dir);
    % get files in this directory
    imgs = strvcat(imgs, spm_get('files', dirn, nw_filter));
  end
  spm_write_sn(imgs,matname,defs.write);
end


