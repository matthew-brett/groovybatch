function [c_ons, c_dur, c_names] = fiac_proc_ref(type, P, TR)
% returns onsets, durations, condition names for fiac experiments

if nargin < 1
  error('Need experiment type');
end
if nargin < 1
  P = spm_get(1, 'fiac*_fonc*.txt');
end
if nargin < 2
  TR = 2.500;
end

c_names = {'SSt_SSp', 'SSt_DSp' 'DSt_SSp', 'DSt_DSp', 'FSt'};
c_ons = {};
c_dur = {};

[ONS COND] = deal(1,2);

ons_data = spm_load(P);
conds = ons_data(:, COND);
ons  = ons_data(:, ONS);

if strcmp(type, 'block')
  % Detect first sentence in block
  cond_change = diff(conds)~=0;
else % event
  % Just pick up first event in session
  cond_change = zeros(length(conds)-1, 1);
end
cond_change = logical([1; cond_change]);
conds(cond_change) = 5;

avg_sentence_dur = 2.277 / TR;

for cno = 1:5
  my_rows = find(conds == cno);
  c_ons{cno} = ons(my_rows) / TR;
  c_dur{cno} = ones(length(my_rows),1) * avg_sentence_dur;
end  

