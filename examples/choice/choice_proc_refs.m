function [c_ons, c_dur, c_parameters] = choice_proc_refs(P, nconds, TR);

% load file
if nargin < 1
  P = spm_get(Inf, '.ref');
end
nrows = size(P,1);
c_ons = cell(nrows, nconds);
c_dur = c_ons;
c_parameters = c_ons;
[COND ONSET DURATION PARAMETER] =deal(1, 2, 3, 4);

for f= 1:size(P,1)
  a = spm_load(deblank(P(f,:)));
  for c = 1:nconds
    my_rows = (a(:,COND) == c);
    c_ons{f, c} = (a(my_rows, ONSET) / 1000) / TR;
    c_dur{f, c} = (a(my_rows, DURATION) / 1000) / TR;
    if size(my_rows, 2) == PARAMETER
      c_parameters{f, c} = a(my_rows, PARAMETER);
    end
  end
end
