function groovy_combine_sessions(glob_ps, sub_ps)
% metabatch file to combine session regressors

% store path
pwd_orig = pwd;
  
for sb = 1:length(sub_ps)
  this_sub = sub_ps(sb);
  % goto SPM results directory
  sub_dir = fullfile(glob_ps.fdata_root, this_sub.dir);
  ana_sdir = fullfile(sub_dir, glob_ps.stats.ana_sdir);
  if ~(exist(ana_sdir))
    error('Have you set up the single subject model?')
  end
  cd(ana_sdir);
  % Get, convert, save design
  load('SPM.mat');
  nSPM = combine_sessions(SPM);
  nSPM.swd = pwd;
  save notcombinedSPM SPM
  savestruct('combinedSPM.mat', struct('SPM', nSPM));
end
cd(pwd_orig);
return

