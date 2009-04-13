function groovy_norm_calc(glob_ps, sub_ps)
% Normalization metabatch; calculates parameters and writes structural only 

defs = glob_ps.defaults.normalise;
for s = 1:length(sub_ps) % for each subject 
  my_sub = sub_ps(s);
  
  % Get normalization source image
  subj(s).P = my_sub.norm.source;
  
  % Make the default normalization parameters file name
  subj(s).matname = [spm_str_manip(subj(s).P,'sd') '_sn.mat'];
  
  % Set the object mask for subject
  subj(s).objmask = my_sub.norm.obj_mask;
  
  % set orientation (.mat file) just in case
  if ~isempty(subj(s).objmask) 
    spm_get_space(subj(s).objmask, ...
		  spm_get_space(subj(s).P));
  end
  
  % Get the images we are going to reslice
  % Because we are going reslice later 
  % We don't reslice anything except the image to be normalized
  subj(s).PP = subj(s).P;      
  
  % call the SPM normalize function to do the work
  spm_normalise(glob_ps.norm.template_images, ...
		subj(s).P, subj(s).matname,...
		defs.estimate.weight, subj(s).objmask, ...
		defs.estimate);
  
  % Do the reslicing
  spm_write_sn(subj(s).PP,subj(s).matname,defs.write);

end











