function groovy_skull_strip(glob_ps, sub_ps)
% batch file to skull strip structural
  
for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  PF = this_sub.segment.source;
  [pn fn ext] = fileparts(PF);
  P = PF;
  for p = 1:3
    P = strvcat(P, fullfile(pn, [fn '_seg' num2str(p) ext]));
  end
  Vi = spm_vol(P);
  Vo = Vi(1);
  Vo.fname = fullfile(pn, [glob_ps.prefixes.skull_stripped fn ext]);
  spm_imcalc(Vi, Vo, 'i1 .* ((i2+i3+i4)>0.5)');
end
