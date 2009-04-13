function groovy_randeff_dir(glob_ps, sub_ps, con_name, rfx_dir)
% Copy contrast images for RFX analysis into their own directory
% FORMAT groovy_randeff_dir(glob_ps, sub_ps, con_name, rfx_dir)
% 
% Last two parameters are optional:
% con_name   - contrast name string, or number, or empty for all contrasts
% rfx_dir    - directory to put contrast images into (empty for default)
%
% If last two not specified, each contrast in the model gets its own
% RFX directory
%
% Matthew Brett 2004/5/19
  
if nargin < 3
  con_name = '';
end
if nargin < 5
  rfx_dir = '';
end
pwd_orig = pwd;
root_dir = glob_ps.fdata_root;

% Disentangle con names 
if ~isempty(con_name) & ischar(con_name)
  con_name = cellstr(con_name);
end
empty_rfx = isempty(rfx_dir);
if ~empty_rfx & ischar(rfx_dir)
  rfx_dir = cellstr(rfx_dir);
end
if ~isempty(con_name) 
  if length(con_name) ~= length(rfx_dir) & ~empty_rfx
    error('Need same number of contrast names and rfx dirs');
  end
else
  if ~isempty(rfx_dir)
    error('Rfx dir should be empty if con_names is empty');
  end
end

% Check if contrasts exists, fill empty con_name, set rfx_dir

% get SPM results directory for first subject
ana_dir = fullfile(root_dir, ...
		   sub_ps(1).dir, ...
		   glob_ps.stats.ana_sdir);

% load SPM model; give "SPM" structure
load(fullfile(ana_dir, 'SPM.mat'));

% Fill empty con_name
if isempty(con_name)
  con_name = 1:length(SPM.xCon);
end

n_cons = length(con_name);
all_c_ns = {SPM.xCon(:).name};

% fill numeric con_name with names
if isnumeric(con_name)
  for c = 1:n_cons
    fname{c} = sprintf('con_%.4d', con_name(c));
  end
  con_name = all_c_ns(con_name);
else
  % Must be cell array - get indices
  for c = 1:n_cons
    xc = SPM.xCon(strmatch(con_name{c}, all_c_ns, 'exact'));
    if isempty(xc)
      error(['Cannot find contrast ' con_name{c}]);
    end
    if ischar(xc.Vcon) % SPM99
      fname{c} = xc.Vcon;
    else               % SPM2
      fname{c} = xc.Vcon.fname;
    end
    fname{c}(end-3:end) = [];
  end
end

% fill rfx_dir

for cn = 1:n_cons
  % filename compatible version of contrast name / no
  fname_con_name{cn} = sf_str2fname(con_name{cn});

  % Fill empty rfx_dir
  if empty_rfx
    rfx_dir{cn} = ['rfx_' fname_con_name{cn}];
  end
  
  % Make rfx_dir absolute path,create if needed
  rfx_dir_full{cn} = fullfile(root_dir, rfx_dir{cn});
  if ~exist(rfx_dir_full{cn}, 'dir')
    mkdir(root_dir, rfx_dir{cn});
  end
end

% Loop across subjects, contrasts, copying images
for sb = 1:length(sub_ps)
  this_sub = sub_ps(sb);

  % get SPM results directory for subject
  ana_dir = fullfile(root_dir, ...
		     this_sub.dir, ...
		     glob_ps.stats.ana_sdir);
  
  for cn = 1:n_cons
    for ext = {'img', 'hdr', 'mat'}
      to_cp = sprintf('%s%c%s.%s', ...
		      ana_dir, ...
		      filesep, ...
		      fname{cn}, ...
		      ext{1});

      if exist(to_cp, 'file')
	cp_to = sprintf('%s%c%s_%s.%s', ...
			rfx_dir_full{cn}, ...
			filesep, ...
			this_sub.dir, ...
			fname_con_name{cn}, ...
			ext{1});
	cmd = sprintf('cp %s %s', to_cp, cp_to);
	fprintf('Copying to %s\n', cp_to);
	unix(cmd);
      end
    end
  end % contrast loop
end % subject loop
return


% Subfunctions

function fname = sf_str2fname(str)

% forbidden chars in file name
badchars = unique([filesep '/\ :;.''"~*?<>|&']);

tmp = find(ismember(str, badchars));   
if ~isempty(tmp)
  str(tmp) = '_';
  dt = diff(tmp);
  if ~isempty(dt)
    str(tmp(dt==1))=[];
  end
end
fname = str;
return
