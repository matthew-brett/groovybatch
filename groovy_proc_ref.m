function [cond_list, c_ons, c_dur, c_parameters] = groovy_proc_ref(filename, sep);
% Load onsets and durations, maybe parametric modulators for session
% FORMAT [conds, c_ons, c_dur, c_parameters] = groovy_proc_ref(filename, sep);
%
% Inputs
% filename - file name containing condition information for subject
% sep      - field separator character (see groovy_load)
%
% Outputs
% conds    - discovered condition names for each row of file
%            (nconds below is length(uniqe(conds)))
% c_ons    - cell array 1 x nconds of onsets for each condition
% c_dur    - cell array 1 x nconds of durations for each condition
% c_parameters - (optional) 
%            cell array 1 x nconds of parametric modulator
%            for each condition
%
% Expects text files with three or four columns, readable by 
% groovy_load.
% Format for each line is condition name, onset, duration, [parameter]
% Values should be separated by char in 'sep' (see groovy_load)
% Onsets and durations are in TRs
% e.g  one line might be 'cond1    34.4   2.3'
% Optionally parametric modulator value as fourth column
  
% load file from GUI if not passed
if nargin < 1
  filename = spm_select(1, '^.*\.ref$', 'Select reference file');
end
if nargin < 2
  sep = [];
end
[conds values] = groovy_load(filename, sep);
[ONSET DURATION PARAMETER] = deal(1, 2, 3);
cond_list = sort(unique(conds));
nconds = length(cond_list);
nrows = 1;
c_ons = cell(nrows, nconds);
c_dur = c_ons;
c_parameters = c_ons;
has_parameters = size(values, 2) == PARAMETER;

for cond_no = 1:nconds
  my_rows = strcmp(conds, cond_list{cond_no});
  c_ons{cond_no} = values(my_rows, ONSET);
  c_dur{cond_no} = values(my_rows, DURATION);
  if has_parameters
    params = values(my_rows, PARAMETER);
    nan_params = isnan(params);
    if any(nan_params)
      if ~all(nan_params)
	msg = sprintf(...
	    'Mixed missing / present params for condition %s, discarding', ...
	    cond_list{cond_no});
	warning(msg);
      end
    else
      c_parameters{cond_no} = params;
    end
  end
end
return
