function groovy_segment(glob_ps, sub_ps)
% batch file to segment structural
  
defs = glob_ps.defaults.segment;

for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  PG = this_sub.segment.template;
  PF = this_sub.segment.source;
  spm_segment(PF, PG, defs);
  
end
