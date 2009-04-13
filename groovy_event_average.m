function [roi_tcs, roi_tcs_names] = groovy_event_average(glob_ps, sub_ps)
% batch file to collect FIR event averages for all events in model
% FORMAT [roi_tcs, roi_tcs_names] = groovy_event_average(glob_ps, sub_ps)
% 
% Returns
% roi_tcs       - ROI fir time courses
% roi_tcs_names - cell array of event names
% 
% Relies on 
%  glob_ps.stats.roi_names
% 
% Warnings
% The events will be in alphabetical order
% Needs marsbar on the path
% 
% $Id: groovy_event_average.m,v 1.2 2005/12/30 11:52:35 matthewbrett Exp $

% Try and start marsbar
marsbar('on')
  
% store path
pwd_orig = pwd;

roi_names = glob_ps.stats.roi_names;
for r = 1:length(roi_names)
  R{r} = maroi(roi_names{r});
end
fir_length = 24;
opts = struct('single', 1, 'percent', 1);

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
  bin_size = tr(D);
  bin_no = fir_length / bin_size;
  
  for r = 1:length(R)
    Y = get_marsy(R{r}, D, 'mean');
    E = estimate(D, Y);
    fir_tc = [];
    for e_t = 1:length(ets)
      fir_tc(:, e_t) = event_fitted_fir(E,...
					ets(e_t).e_spec, ...
					bin_size, ...
					bin_no, ...
					opts);
      
    end
    roi_tcs{s, r} = fir_tc;
    roi_tcs_names = {ets(:).name};
  end
  
  
end

% back to initial directory
cd(pwd_orig);
