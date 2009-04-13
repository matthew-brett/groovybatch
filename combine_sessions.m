function nSPM = combine_sessions(SPM)
% Makes event model combined across sessions
% FORMAT nSPM = combine_sessions(SPM)
%
% Inputs
% SPM     - spm FMRI design structure (with event and data definitions)
%
% Output
% nSPM    - modified design with event regressors for sessions combined
%
% Typical usage
% >> load SPM;
% >> mkdir combined
% >> cd combined
% >> SPM = combine_sessions(SPM);
% >> SPM = spm_spm(SPM); % estimate
%
% Event columns are concatenated across sessions
% Non-event regressors (such as movement parameters) are left as session
% specific.  SPM has row-specific filters and autocorrelation estimation,
% so the sessions will continue to have correct (separate) filtering and
% autocorrelation estimates after combination.
%
% Note that combining across sessions has the effect of leaving between
% session variance in the variance estimates, maybe a good thing, maybe
% not. 
%
% Tested for SPM5, might work for SPM2
%
% Needs marsbar on your matlab path

try
  marsbar on
catch
  error('You need marsbar on the path')
end
 
D = mardo(SPM);
if ~is_fmri(D)
  error('Marsbar does not recognize this as an FMRI design');
end
if ~has_images(D)
  error('Need design with data');
end
br = block_rows(D);
X = design_matrix(D);
XM = size(X, 2);
ev_col_list = sf_session_event_cols(D);
old_names = SPM.xX.name;
hrfXs = {};
hrf_names = {};
for ss = 1:length(br)
  ev_cols = ev_col_list{ss};
  hrfXs{ss} = X(br{ss}, ev_cols);
  hrf_names{ss} = old_names(ev_cols);
end
[newXhrf, hrf_names] = sf_stack_hrfs(hrfXs, hrf_names);
other_cols = 1:XM;
other_cols(ismember(other_cols, [ev_col_list{:}])) = [];
newX = [newXhrf X(:,other_cols)];
names = [hrf_names old_names(other_cols)];
% Now make the new design
X_col_diff = XM - size(newX, 2);
xX.X = newX;
xX.name = names;
% Assuming block and global effects are after event columns
xX.iB = SPM.xX.iB - X_col_diff;
xX.iG = SPM.xX.iG - X_col_diff;
xX.K = SPM.xX.K;
nSPM = struct(...
    'SPMid', SPM.SPMid,...
    'xY', SPM.xY,...
    'xX', xX,...
    'xVi', SPM.xVi,...
    'xM', SPM.xM);
return

function col_list = sf_session_event_cols(D)
% Return the event modeling columns for each session
all_es = event_specs(D);
session_numbers = all_es(1,:);
if any(diff([0 session_numbers])>1)
  error('Out of order sessions');
end
for ss = unique(session_numbers)
  es = all_es(:,all_es(1,:)==ss);
  cols = [];
  for esno = 1:size(es,2)
    cols = [cols event_cols(D, es(:,esno))];
  end
  if ~all(cols == sort(unique(cols)))
    error('Non-unique or out of order columns');
  end
  col_list{ss} = cols;
end
return

function [s_hrf, c_names] = sf_stack_hrfs(hrfXs, hrf_names)
% Strip session identifiers, find unique
[x_names, c_names] = sf_strip_session(hrf_names);
% Across sessions, fill any present regressors, by name
s_hrf = [];
n_conds = length(c_names);
for ss = 1:length(hrfXs)
  hrfX = hrfXs{ss};
  ss_names = x_names{ss};
  nX = zeros(size(hrfX, 1), n_conds);
  for cno = 1:n_conds
    matching_i = ismember(ss_names, c_names{cno});
    if any(matching_i)
      nX(:, cno) = hrfX(:, matching_i);
    end
  end
  s_hrf = [s_hrf; nX];
end
return

function [new_names, unique_names] = sf_strip_session(names)
% Strip session identifiers from event column names
new_names = {};
all_names = {};
for ss = 1:length(names)
  ss_id = sprintf('Sn(%d) ', ss);
  first_i = length(ss_id)+1;
  ss_names = names{ss};
  nss_names = {};
  for nno = 1:length(ss_names)
    nss_names{nno} = ss_names{nno}(first_i:end);
  end
  new_names{ss} = nss_names;
  all_names = [all_names nss_names];
end
unique_names = sort(unique(all_names));
return
