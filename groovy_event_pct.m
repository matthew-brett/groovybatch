function [roi_pct, roi_pct_names] = groovy_event_pct(glob_ps, sub_ps)
% batch file to collect percent signal changes for all events in model
% FORMAT [roi_pct, roi_pct_names] = groovy_event_pct(glob_ps, sub_ps)
% 
% Returns
% roi_pct       - percent signal change for each event type in model
% roi_pct_names - cell array of event names
% 
% Relies on 
%  glob_ps.stats.roi_names
% 
% Warnings
% The events will be in alphabetical order
% Needs marsbar on the path
% 
% $Id: groovy_event_pct.m,v 1.2 2005/12/30 11:52:35 matthewbrett Exp $

% Try and start marsbar
marsbar('on')
  
% store path
pwd_orig = pwd;

roi_names = glob_ps.stats.roi_names;
dur = glob_ps.stats.roi_dur;
for r = 1:length(roi_names)
  R{r} = maroi(roi_names{r});
end

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
  
  D = mardo(SPM);
  ets = event_types_named(D);
  
  for r = 1:length(R)
    Y = get_marsy(R{r}, D, 'mean');
    E = estimate(D, Y);
    pct = [];
    for e_t = 1:length(ets)
      pct(e_t) = event_signal(E,...
			      ets(e_t).e_spec, ...
			      dur);
      
    end
    roi_pct{s, r} = pct;
    roi_pct_names = {ets(:).name};
  end
  
  
end

% back to initial directory
cd(pwd_orig);
