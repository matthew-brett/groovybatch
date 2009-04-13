function [g_params, s_params] = fiac_top_groove(subjects, g_params, flags)
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

% flags are specific to this (FIAC) experiment
def_flags = struct('experiment_types', 'all', ...
		   'fm_process', 'real_imag');
  
all_subjects = {...
    'fiac1',...
    'fiac2',...
    'fiac3',...
    'fiac4',...
    'fiac6',...
    'fiac7',...
    'fiac9',...
    'fiac10',...
    'fiac12',...
    'fiac13',...
    'fiac14',...
    'fiac15'...
	       };

% Subject exclusion:
% fiac0 - no fieldmap
% fiac11 - no fieldmap
% fiac5 excluded - no anatomical
% fiac8 excluded - massive spikes

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
  g_params.fdata_root = '/home/imagers/FIAC';
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

% Check the default orientation is as the FIAC images expect
if ~spm_flip_analyze_images
  error('Please check SPM defaults _do_ flip analyze images for FIAC data');
end

% prefixes for processed images (prepended to raw image filter below, in
% the scripts)
g_params.prefixes = struct( ...
    'slice'         , '', ...    % prefix for images _to_ do slice timing on
    'realign'       , 'a', ...   % prefix to select imgs _to_ realign 
    'coreg'         , 'ua', ...  % prefix for images _to_ coregister
    'skull_stripped', 'ss_', ... % prefix for skull stripped images
    'norm_write'    , 'ua', ...
    'smooth'        , 'wua', ...
    'stats'         , 'swua');

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
% Unwarp 
% -------------------------------------------

% unwarp parameters = see spm_realign_ui for definitions
g_params.fieldmap.foe = [4 5];  % first order effects
g_params.fieldmap.soe = [];     % second order effects
g_params.fieldmap.do_realign = 0; % whether to do realignment during unwarp

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
g_params.stats.FWHM = 5; % mm

% -------------------------------------------
% Statistics
% -------------------------------------------

% subdirectory name for analysis
g_params.stats.ana_sdir = 'spm2_ana_5mm';

% Explicit mask - usually missing or empty
% In our case, set to mask of gray+white+csf
mskdir = 'masks';
mskimg = 'anything.img';
mskfull = fullfile(g_params.fdata_root, mskdir);
if ~exist(mskfull, 'dir')
  mkdir(g_params.fdata_root, mskdir);
end
g_params.stats.explicit_mask = fiac_make_mask(...
    fullfile(mskfull, mskimg), g_params.defaults);

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
% Phiwave
% -------------------------------------------
if ~isempty(which('phiwave'))
  phiwave('on');
  scales = 4;

  % Set some options for estimation (see @phido/estimate.m)
  g_params.phiwave.estimate =  struct(...
      'wavelet',  phiw_lemarie(2), ...
      'scales',   scales, ...
      'write_res', 1, ...
      'wtprefix', 'wv_');
  
  % Set some denoising options (see @phido/get_wdimg.m)
  g_params.phiwave.denoise = struct(...
      'thcalc', 'stein', ...
      'thapp', 'linear', ...
      'levels', [ones(1, scales) -Inf], ...
      'write_err', 1);
  
  g_params.phiwave.contrasts = struct(...
      'suffix',  '_stein');
end


% -------------------------------------------
% Parameters used in subject / session loop
% -------------------------------------------

% in seconds
TR = 2.5; 

% Whether to store structural stuff on parameter or fdata root
my_anat_root = g_params.parameter_root;

% Session directories, in fact same for each subject
my_sesses = {...
    'fonc1',...
    'fonc2',...
    'fonc3',...
    'fonc4'...
	    };

% Transformations to apply to all images per subject before we start
% Cam be different for each subject.  Can be vector (translations,
% rotations, zooms).
% Can also be a 4x4 transformation, in which case the orientation for
% each image is set to be this 4x4 matrix.
trans = 0;

% Time to acquire one slice of data - same for each subject here
slice_time = 0.083;

% Slice acquisition order - same for each subject here
acq_order = [1:2:30 2:2:30];

% Fieldmap defaults, from subroutine below, with comments
% Can be different for each subject; some setting are the same, though
pm_def = sf_pm_settings;
pm_def.match_vdm = 0;

% covariate names - in fact same for each subject and session
% For none, need empty cell array "{}"
cov_names = {};

% Number of columns we will use for movement parameters (for contrast
% construction in subject loop)
n_m_params = length(g_params.stats.movement_params) * 6;

% Contrast names for both experiments
con_names = {'FSt',...
	     'Dst-SSt', ...
	     'DSp-SSp',...
	     'Interaction'...
	    };


% Here one row per contrast, one column per condition in each (and every)
% session
con_mat = [ 0  0  0  0  1; ...
	   -1 -1  1  1  0; ...
	   -1  1 -1  1  0; ...
	    1 -1 -1  1  0 ...
	  ];

if ~isempty(which('phiwave'))
  g_params.phiwave.contrasts.numbers = 2:size(con_mat, 1)+1;
end

% Decode session selecting flags (specific to FIAC experiment)
[get_ev_f get_blk_f] = sf_get_expt_types(flags);

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
  
  % No first run for this subject
  if strcmp(sub_str, 'fiac7')
    sub_sesses(1) = [];
  end
  
  % Subject was asleep for 3rd session 
  if strcmp(sub_str, 'fiac10')
    sub_sesses(3) = [];
  end
  
  % Remove sessions with big variance spikes
  if strcmp(sub_str, 'fiac11')
    sub_sesses(3) = [];
  end
  if strcmp(sub_str, 'fiac12')
    sub_sesses(1) = [];
  end
  
  n_sesses = length(sub_sesses);
  
  % Filter for raw image names.
  sub_struct.raw_filter = [sub_str '_fonc*.img'];
  
  % Any transformations to apply for this subject
  % (set to 0 if none required)
  sub_struct.trans = 0;

  % TR for each subject.  Sometimes it's different for each subject
  % but in this case it's the same
  sub_struct.TR = TR;

  % Slice timing stuff
  sub_struct.slice.time = slice_time;
  sub_struct.slice.acq_order = acq_order;
  sub_struct.slice.ref_slice = 1;
  
  % Segment stuff
  anat_name = sprintf('%s_anat1', sub_str);
  anat_dir = fullfile(my_anat_root, sub_str, 'anat1');
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
  cond_dir = sub_dir_p;

  % fieldmap directory, and file ends
  fm_dir = fullfile(my_anat_root, sub_str, 'fieldmap2');
  fm_ends = {'0_real', '0_imag', '2_real', '2_imag'};
  
  switch flags.fm_process
    case 'real_imag'
     % This to do all fieldmap processing, from real/imag image pairs
     P = [];
     for i = 1:length(fm_ends)
       fname = sprintf('%s_fieldmap2_acq%s.img', sub_str, fm_ends{i});
       P = strvcat(P, fullfile(fm_dir, fname));
     end
   case 'vdm_only'
    % This, just to create the VDM map from a pre-existing fieldmap
    fname1 = sprintf('%s_fieldmap2_acq%s.img', sub_str, fm_ends{1});
    P = fullfile(fm_dir, ['fpm_' fname1]);
    P = strvcat(P, fullfile(fm_dir, ['mag_' fname1]));
   otherwise
    error(['Odd fieldmap processing: ' flags.fm_process]);
  end
  
  % We need to count up number of valid sessions to allow check if this
  % is a block or event related design
  ss_ctr = 1;
  
  % -------------------------------------------
  % Session loop
  % -------------------------------------------
  
  for ss = 1:n_sesses
    % Session structure, within subject structure
    ss_struct = [];
    
    % First check if we are processing this session, by checking cond
    % file name.  This is specific to the way the batching works for the
    % two-experiments-per-subject thing of FIAC
    sub_no = sscanf(sub_str, 'fiac%d');
    cond_file = spm_get('files', ...
			cond_dir, ...
			sprintf('subj%d_*_%s.txt', sub_no, sub_sesses{ss}));

    if isempty(cond_file), 
      error('No condition files for sub %s, session %s', ...
	     sub_str, sub_sesses{ss})
    end
    if ~isempty(findstr(cond_file, '_bloc_')) % this is a block session
      if ~get_blk_f, continue,  end  % abort processing this session
      [ss_struct.ons ss_struct.dur ss_struct.cond_names] = ...
	  fiac_proc_ref('block', cond_file, sub_struct.TR);
    elseif ~isempty(findstr(cond_file, '_evt_')) % this an event session
      if ~get_ev_f, continue, end % abort processing this session
      [ss_struct.ons ss_struct.dur ss_struct.cond_names] = ...
	  fiac_proc_ref('event', cond_file, sub_struct.TR);
    else
      error(['Odd stimulus file ', cond_file]);
    end
    
    % Field map stuff
    ss_struct.fieldmap = pm_def;
    ss_struct.fieldmap.imgs  = P;
    fname1 = sprintf('vdm_%s_fieldmap2_acq%s.img', sub_str, fm_ends{1});
    ss_struct.fieldmap.vdm_img = fullfile(fm_dir, fname1);

    % Fill session structure
    ss_struct.dir = sub_sesses{ss}; % here the directory names are all the same
    ss_struct.covs = [];
    ss_struct.cov_names = cov_names;
    
    % Put into subject structure
    sub_struct.sesses(ss_ctr) = ss_struct;
    ss_ctr = ss_ctr + 1;
  end
  n_sesses = ss_ctr - 1;
  
  % For the coregistration steps, to know the name of the eventual mean
  % image, we need to know the first image in the series, which will
  % probably be called image_0005.img or something.  The step below could
  % also be done with spm_get('files', and the raw_filter field.
  fname = sprintf('%s_%s_0005.img', sub_str, sub_struct.sesses(1).dir);
  first_img = fullfile(g_params.fdata_root, ...
		       sub_str, ...
		       sub_struct.sesses(1).dir, ...
		       fname);
  
  % Fieldmap target (EPI) image
  sub_struct.fieldmap.target = sf_img_prefix(...
      first_img, ...
      ['mean' g_params.prefixes.realign]);
  
  % Coreg target image
  sub_struct.coreg.target = sf_img_prefix(...
      first_img, ...
      ['mean' g_params.prefixes.coreg]);
  
  % And to-be-coregistered image
  sub_struct.coreg.object = fullfile(anat_dir, ...
					 [anat_name '.img']);
  
  % And images to take along with coregistration
  % dirs field here has to be cell array to cope with many dirs, for
  % sessions (if 'other' are EPI images for example)
  sub_struct.coreg.other_dirs = {anat_dir};
  sub_struct.coreg.other_filter = [anat_name '_seg*.img'];
  
  % Put contrasts
  n_cons = size(con_mat, 1);
  ss_con = [con_mat, zeros(n_cons, n_m_params)];
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

% --------------------------------------------------
% Subfunctions specific to FIAC experiment
% --------------------------------------------------

function [get_ev_f, get_blk_f] = sf_get_expt_types(flags)
% subfunction to process flags to select experiment types
  
get_ev_f = 1;
get_blk_f = 1;
if ~isfield(flags, 'experiment_types'), return, end
f_t = flags.experiment_types;
if isempty(f_t), return, end
if ischar(f_t), f_t = {f_t}; end
if ismember('all', f_t), return, end
get_blk_f = ismember('block', f_t);
get_ev_f = ismember('event', f_t);
return

function fiac_pm = sf_pm_settings
% Sets the values for the FieldMap toolbox for FIAC dataset

% Defaults for creating field map. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
fiac_pm.INPUT_DATA_FORMAT = 'RI';      % 'RI' = load two real and 
                                      % imaginary image pairs
                                      % 'PM' = load one or two
                                      % phase and magnitude image
                                      % pairs.
fiac_pm.SHORT_ECHO_TIME = 4.2;        % Short echo time in ms
fiac_pm.LONG_ECHO_TIME = 13.304;      % Long echo time in ms

% Defaults for unwrapping options. (See pm_make_fieldmap.m and 
%                                   FieldMap.man for more info.)
%=======================================================================
fiac_pm.UNWRAPPING_METHOD = 'Mark3D';  % Unwrapping options are:
                                      % 'Huttonish', 'Mark3D' or 'Mark2D'
fiac_pm.FWHM = 10;                     % FWHM of Gaussian filter used to 
                                      % implement weighted smoothing of
                                      % unwrapped maps.
fiac_pm.PAD = 10;                       % Size of padding kernel if required.
fiac_pm.WS = 0;                        % Weighted or normal smoothing.

% Defaults for converting field map to voxel displacement map.
%=======================================================================
fiac_pm.EPI_BASED_FIELDMAPS = 0;          % EPI=1, other=0.
fiac_pm.K_SPACE_TRAVERSAL_BLIP_DIR = -1;  % +ve k-space = 1, -ve = -1.
fiac_pm.TOTAL_EPI_READOUT_TIME = 33.1776; % Sonata EPI RO time (500E-6*64)

% Defaults for Unwarping.
%=======================================================================
fiac_pm.DO_JACOBIAN_MODULATION = 1;   % Do jacobian modulation to adjust 
                                      % for compression or stretching
                                      % No = 0, Yes = 1

% And for graphics
fiac_pm.SHOW_GRAPHICS = 0;            % Whether to show graphics 
return

