function groovy_smooth(glob_ps, sub_ps)
% metabatch spm2 smooth script

imgs = '';
for sb = 1:length(sub_ps)
  my_sub = sub_ps(sb);
  subjroot = fullfile(glob_ps.fdata_root, my_sub.dir);
  sm_filter = [glob_ps.prefixes.smooth my_sub.raw_filter];
  for ss = 1:length(my_sub.sesses)
    dirn = fullfile(subjroot,my_sub.sesses(ss).dir);
    % get files in this directory
    imgs = strvcat(imgs, spm_get('files', dirn, sm_filter));
  end
end
  
% and this is just spm_smooth_ui pasted/edited
s     = glob_ps.stats.FWHM;
P     = imgs;
n     = size(P,1);

% implement the convolution
%---------------------------------------------------------------------------
for i = 1:n
  Q = deblank(P(i,:));
  [pth,nm,xt,vr] = fileparts(deblank(Q));
  U = fullfile(pth,['s' nm xt vr]);
  spm_smooth(Q,U,s);
end




