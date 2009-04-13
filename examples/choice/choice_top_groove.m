function [g_params, s_params] = choice_top_groove(subjects, g_params, flags)
% sets up various things specific to this analysis
% 
% FORMAT [g_params, s_params] = this_function(subjects, g_params, flags)
%   
% Input
% subjects       - subject directory name identifying subject OR
%                  cell array of subject directory names OR
%                  vector of subject numbers
%                  If empty, all subjects' data returned
% g_params  - structure overriding default things like
%                  'fdata_root' root directory for functional data
%                  'parameter_root' root directory for parameters
% flags     - any interesting flags to change behaviour of the script
%
% Outputs
% g_params  - parameters shared across all subjects' analyses
% s_params - struct array of parameters per subject
%
% The fields in g_params contain all the stuff that is common for all
% subjects - like the session names.
%
% The s_params struct array contains one struct per subject, with one
% field ('dir') giving the subject directory name - e.g s_params(1) =
% struct('dir', '00AH', 'some_things', 1, 'other_things', 2) and so on.
% This way we dont necessarily have to keep track of things like subject
% _numbers_ we can reference the subjects details with a string - e.g. to
% get the TR for a particular subject:
% 
% [gp sp] = this_function('00AH');
% TR = sp.TR
% 
% OR
%  
% [gp sp] = this_function;
% sub_no = strmatch('00AH', {sp(:).dir}, 'exact');
% TR = sp(sub_no).TR
% 
% First written by MB, 30 March 2005

% flags are specific to this (choice) experiment
def_flags = struct('svd', 0);
  
all_subjects = {...
    '03FR',...
    '04AL',...
    '05AM',...
    '06AW',...
    '07CS',...
    '09NM',...
    '10PL',...
    '11KS',...
    '12KK',...
    '13AM',...
    '14AD',...
    '15BD',...
    '16NK',...
    '18RL'...
	       };

% Subject excluded - massive artefacts '17HI',...

if nargin < 1
  subjects = [];
end
if nargin < 2
  g_params = [];
end
if nargin < 3
  flags = [];
end
flags = groovy_struct('ffillsplit', def_flags, flags);

if isempty(subjects)
  subjects = all_subjects;
elseif isnumeric(subjects)
  subjects = all_subjects(subjects);
elseif ischar(subjects)
  subjects = cellstr(subjects); 
end
nsubs = length(subjects);

% Where the subjects' data directories are stored
if ~isfield(g_params, 'fdata_root')
  g_params.fdata_root = '/home/imagers/choice';
end

% Where the directory tree containing parameters is stored 
% This could be where you have stored your batch files, some already
% calculated normalization parameters, reference functions and so on.
% parameter_root could be the same as fdata_root above.
% See the TWiki
% http://dynevor.hopto.org/twiki/bin/view/IvryImaging/HowtoHummer
% for more explanation.  
if ~isfield(g_params, 'parameter_root')
  g_params.parameter_root = g_params.fdata_root;
end

% Get SPM defaults, so we can modify them 
spm_defaults;
g_params.defaults = defaults;

% prefixes for processed images (prepended to raw image filter below, in
% the scripts)
g_params.prefixes = struct( ...
    'slice'         , '', ...    % prefix for images _to_ do slice timing on
    'realign'       , '', ...   % prefix to select imgs _to_ realign 
    'coreg'         , '', ...  % prefix for images _to_ coregister
    'skull_stripped', 'ss_', ... % prefix for skull stripped images
    'norm_write'    , '', ...
    'smooth'        , 'w', ...
    'stats'         , 'sw');

% -------------------------------------------
% Realignment
% -------------------------------------------
% what to output from the realignment
%                  0   - don't create any resliced images.
%                        Useful if you only want a mean resliced image.
%                  1   - don't reslice the first image.
%                        The first image is not actually moved, may not be
%                        necessary to resample it.
%                  2   - reslice all the images.
g_params.realign.which = 2;  % needed for unwarp output

% Whether or not to write mean image
g_params.realign.write_mean = 1;

% Default parameter changes to apply to realignment
g_params.defaults.realign.write.interp = 1;  % trilinear interpolation

% -------------------------------------------
% Coregistration
% -------------------------------------------

% which method to use for coregistration ('spm','flirt','flirt_9')
g_params.coreg.method = 'spm';
				   
% -------------------------------------------
% Normalization
% -------------------------------------------

% Turn off template weighting
g_params.norm.estimate.weight = '';

% list of templates for normalization
% Here, a smoothed grey matter segmentation of MNI brain 
g_params.norm.template_images = ...
    fullfile(spm('Dir'),'apriori', 'gray.mnc');

% -------------------------------------------
% Smoothing
% -------------------------------------------
g_params.stats.FWHM = 8; % mm

% -------------------------------------------
% Statistics
% -------------------------------------------

% subdirectory name for analysis
g_params.stats.ana_sdir = 'spm2_humm_ana';

if flags.svd
  g_params.stats.ana_sdir = [g_params.stats.ana_sdir '_svd'];
end

% Common stuff about experimental design here 

% Event basis function
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
g_params.stats.event_bf.name = 'hrf';
g_params.stats.event_bf.length = 20;  % this isn't used for an hrf model
g_params.stats.event_bf.order = 1;    % this isn't used for an hrf model

% Units for onsets ('scans'|'secs')
g_params.stats.units = 'scans';

% High pass filter in seconds
g_params.stats.high_pass = 120;

% intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w'
g_params.stats.ar_form = 'none';

% What movement parameters to include in the model
% Cell array of one or more of 
% 'moves'   - user Nx^ movement parameters as regressors
% 'mm1'     - use movement parameters lagged by 1
% 'moves_2' - movement parameters squared
% 'mm1_2'   - mm1 squared
% (empty cell array = no movement parameter correction)
g_params.stats.movement_params = {'moves'};

% -------------------------------------------
% Parameters used in subject / session loop
% -------------------------------------------

% in seconds
TR = 1.934; 

% Whether to store structural stuff on parameter or fdata root
my_anat_root = g_params.parameter_root;

% Session directories, in fact same for each subject
my_sesses = {...
    'epi1',...
    'epi2',...
    'epi3',...
    'epi4',...
    'epi5',...
    'epi6',...
	    };

% Transformations to apply to all images per subject before we start
% Cam be different for each subject.  Can be vector (translations,
% rotations, zooms).
% Can also be a 4x4 transformation, in which case the orientation for
% each image is set to be this 4x4 matrix.
trans = 0;

% Time to acquire one slice of data - same for each subject here
slice_time = NaN;

% Slice acquisition order - same for each subject here
acq_order = [1:2:30 2:2:30];

% condition stuff - in fact same for each subject and session 
cond_names = {'Direct4', 'Direct2', 'Choice', 'Symbolic'};

% covariate names - in fact same for each subject and session
% For none, need empty cell array "{}"

cov_names = {'eig1',...
	     'eig2',...
	     'eig3',...
	     'eig4',...
	     'eig5'};

% Number of columns we will use for movement parameters (for contrast
% construction in subject loop)
n_m_params = length(g_params.stats.movement_params) * 6;

% Contrasts
con_names = {'All moves', ...
	     'Direct4', ...
	     'Direct2',...
	     'Choice', ...
	     'Symbolic',...
	     'Symbolic - Direct2', ...
	     'Symbolic - Direct4',...
	     'Choice - Direct2', ...
	     'Choice - Direct4'};

    
con_mat = [1 1 1 1; ...
	   1 0 0 0; ...
	   0 1 0 0; ...
	   0 0 1 0; ...
	   0 0 0 1; ...
	   0 -1 0 1; ...
	   -1 0 0 1; ...
	   0 -1 1 0; ...
	   -1 0 1 0];

% -------------------------------------------
% Subject / session loop
% -------------------------------------------

% Now set subject by subject stuff
for sb = 1:nsubs
  clear sub_struct
  sub_str = subjects{sb};
  sub_dir_f = fullfile(g_params.fdata_root, sub_str);
  sub_dir_p = fullfile(g_params.parameter_root, sub_str);
  sub_sesses = my_sesses;
  
  % -------------------------------------------
  % Session exclusion
  % -------------------------------------------
  
  % Subject 05 didn't get a session 6
  if strcmp(sub_str, '05AM')
    sub_sesses(end) = [];
  end
  
  % Subject 09 does not have a session 5
  if strcmp(sub_str, '09NM')
    sub_sesses(5) = [];
  end
  n_sesses = length(sub_sesses);
  
  % Filter for raw image names.
  sub_struct.raw_filter = [sub_str '_epi*.img'];
  
  % Any transformations to apply for this subject
  % (set to 0 if none required)
  sub_struct.trans = trans;

  % TR for each subject.  Sometimes it's different for each subject
  % but in this case it's the same
  sub_struct.TR = TR;

  % Slice timing stuff
  sub_struct.slice.time = slice_time;
  sub_struct.slice.acq_order = acq_order;
  sub_struct.slice.ref_slice = 1;
  
  % Segment stuff
  anat_name = 'mpflash';
  anat_dir = fullfile(my_anat_root, sub_str);
  sub_struct.segment.source = fullfile(anat_dir, ...
				       [anat_name '.img']);
  sub_struct.segment.template = fullfile(spm('dir'), ...
					 'templates', ...
					 'T1.mnc');
  
  % image to normalize for this subject
  sub_struct.norm.source = fullfile(anat_dir, ...
				    [anat_name '_seg1.img']);

  % Other images (in space of structural) to write normalized
  % (epis resliced in another write-normalized pass)
  sub_struct.norm.others = fullfile(anat_dir, ...
				    [anat_name '.img']);
  
  % object mask - usually empty
  sub_struct.norm.obj_mask = '';
  
  % condition information.
  % read from one file for each session
  % Also covariates (here not used)
  cond_dir = fullfile(sub_dir_p, 'refs');

  if flags.svd
    % Load previous SPM, SVD
    svd_file = fullfile(g_params.fdata_root, ...
			sub_str, ...
			'spm2_humm_ana', ...
			'SVD.mat');
    load(svd_file);
    eigs = SVD.u(:,1:5);
    clear SVD
    spm_file = fullfile(g_params.fdata_root, ...
			sub_str, ...
			'spm2_humm_ana', ...
			'SPM.mat');
    load(spm_file);
    row = {SPM.Sess(:).row};
    clear SPM
  end
  
  % -------------------------------------------
  % Session loop
  % -------------------------------------------
  
  for ss = 1:n_sesses
    % Session structure, within subject structure
    ss_struct = [];
    
    % Condition file
    cond_file = spm_get('files', cond_dir, ...
			sprintf('%s_%s.ref', sub_str, sub_sesses{ss}));
    
    if isempty(cond_file), 
      error('No condition files for sub %s, session %s', ...
	     sub_str, sub_sesses{ss})
    end

    % Fill session structure
    ss_struct.dir = sub_sesses{ss}; % here the directory names are all the same
    [ss_struct.ons ss_struct.dur ss_struct.parameters] = ...
	choice_proc_refs(cond_file, length(cond_names), sub_struct.TR);
    if isempty(parameters)
      ss_struct.P = struct('name', 'none');
    else
      ss_struct.P = struct('name', 'other',...
			   'P', parameters(:), ...
			   'h', 1); % h is order, 1 being just linear
    end
    ss_struct.cond_names = cond_names;

    if flags.svd
      ss_struct.covs = eigs(row{ss}, :);
      ss_struct.cov_names = cov_names;
    else
      ss_struct.covs = [];
      ss_struct.cov_names = {};
    end
    
    % Put into subject structure
    sub_struct.sesses(ss) = ss_struct;
  end
  
  % For the coregistration steps, to know the name of the eventual mean
  % image, we need to know the first image in the series, which will
  % probably be called image_0001.img or something.  The step below could
  % also be done with spm_get('files', and the raw_filter field.
  fname = sprintf('%s_%s_0001.img', sub_str, sub_struct.sesses(1).dir);
  first_img = fullfile(g_params.fdata_root, ...
		       sub_str, ...
		       sub_struct.sesses(1).dir, ...
		       fname);
  
  % Coreg target image
  sub_struct.coreg.target = sf_img_prefix(...
      first_img, ...
      ['mean' g_params.prefixes.coreg]);
  
  % And to-be-coregistered image
  sub_struct.coreg.object = fullfile(anat_dir, [anat_name '.img']);
  
  % And images to take along with coregistration
  % dirs field here has to be cell array to cope with many dirs, for
  % sessions (if 'other' are EPI images for example)
  sub_struct.coreg.other_dirs = {anat_dir};
  sub_struct.coreg.other_filter = [anat_name '_seg*.img'];
  
  % Put contrasts
  n_cons = size(con_mat, 1);
  n_covs = n_m_params + (5 * flags.svd);
  ss_con = [con_mat, zeros(n_cons, n_covs)];
  sb_con = [repmat(ss_con, 1, n_sesses) zeros(n_cons, n_sesses)];
  sub_struct.contrasts = struct(...
      'name', con_names',...
      'value', num2cell(sb_con, 2), ...
      'type', 'T');
  
  % Set the subdirectory for this subject.  This is required 
  sub_struct.dir = sub_str;
  
  % Set into returned structure
  s_params(sb) = sub_struct;
  
end


return

% --------------------------------------------------
% Generic subfunctions
% --------------------------------------------------

function fname = sf_img_prefix(fname, prefix)
[pn fn ext] = fileparts(fname);
fname = fullfile(pn, [prefix fn ext]);
return

