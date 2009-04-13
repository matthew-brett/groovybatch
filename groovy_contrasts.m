function groovy_contrasts(glob_ps, sub_ps)
% batch file to run contrasts

% store path
pwd_orig = pwd;

for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  % get, goto SPM results directory
  ana_dir = fullfile(glob_ps.fdata_root, ...
		     this_sub.dir, ...
		     glob_ps.stats.ana_sdir);
  cd(ana_dir);
  
  % load SPM model; give "SPM" structure
  disp('Loading SPM.mat');
  load('SPM.mat');
  disp('Done');
  
  % Fix swd, just in case
  SPM.swd = ana_dir;
  
  % Where we will start filling in the contrasts
  xcon_1 = length(SPM.xCon)+1;
  
  % now put contrast into SPM structure
  for cn = 1:length(this_sub.contrasts)
    c = this_sub.contrasts(cn);
    if c.type == 'T'  & any(size(c.value) == 1)
      c.value = c.value(:);
    end
    SPM.xCon(end + 1)= spm_FcUtil('Set',...
				  c.name,...
				  c.type,...
				  'c',...
				  c.value, ...
				  SPM.xX.xKXs);
  end
  
  % Estimate only the contrasts we've added
  spm_contrasts(SPM, xcon_1:length(SPM.xCon));
  
end

% back to initial directory
cd(pwd_orig);
