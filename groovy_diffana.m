function groovy_diffana(glob_ps, sub_ps)
% time series difference analysis metabatch file
  
% load spm defaults; we need to store them, then set the print command
% into the global to change the print output.
spm_defaults;
global defaults;
d = defaults;

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  ps_file = [this_sub.dir '_tsplot.ps'];
  defaults.printstr = [spm_figure('DefPrintCmd'), ps_file];
  s_filter = this_sub.raw_filter;
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_get('files', dirn, s_filter);
 
    % do analysis
    tsdiffana(P);
    
    % print
    spm_figure('print');
  end
end
defaults = d;  
