function groovy_svd(glob_ps, sub_ps)
% metabatch file to run SVD on model
%
% To run SVD, we need to make an argfile to run the MM analysis.  We need
% these:
% [gcsdr cwd csd d fname ImgDir ImgName paramsAnal.temporalFilter paramsAnal.divideByRessd paramsAnal.resContSp] 
% The parameters required are (in order, one line per subject)
% gcsd		= global result directory.
% cwd		= current working directory (dir contening the model)
% csd 		= current saving directory for each space domain.
% d             = contrast number
% fname 	= Eigen image basename.
% ImgDir        = Directory containing images
% ImgName       = Filter to select image in ImgDir
% paramsAnal.temporalFilter = whether to apply analysis temporal filter
% paramsAnal.divideByRessd  = whether to divide each voxel by voxel variance
% paramsAnal.resContSo = whether to project into residual space or not

if nargin < 1
  glob_ps = '';
end
if nargin < 2
  sub_ps = '';
end
if isempty(glob_ps)
  glob_ps = spm_get([0 1], 'SPM.mat', 'Select analysis for SVD');
end
if isempty(glob_ps), return, end

% glob_ps can be an SPM.mat filename
if ischar(glob_ps)
  if ~isempty(sub_ps)
    error('Need missing or empty sub_ps with string glob_ps');
  end
  str = glob_ps;
  clear sub_ps glob_ps;
  str = fileparts(str);
  [str glob_ps.stats.ana_sdir] = fileparts(str);
  [glob_ps.fdata_root sub_ps.dir ] = fileparts(str);
end

% default parameters
svd_params = struct(...
    'contrast_no', 1, ...
    'eigen_name' , 'eigen_', ...
    'temporal_filter', 1, ...
    'normalize', 0, ...
    'project_residual', 0);

% replace any defaults with stuff from glob_ps svd struct
if isfield(glob_ps, 'svd')
  svd_params = groovy_struct('ffillsplit', svd_params, glob_ps.svd);
end

for s = 1:length(sub_ps) % for each subject 
  this_sub = sub_ps(s);
  
  % get, goto SPM results directory
  ana_dir = fullfile(glob_ps.fdata_root, ...
		     this_sub.dir, ...
		     glob_ps.stats.ana_sdir);
  
  argfile = tempname;
  fid = fopen(argfile, 'wt');
  if fid == -1, error(['Cannot open file ' tmp]); end
  fprintf(fid, '%s %s %s %d %s %s %s %d %d %d', ...
	  ana_dir, ...
	  ana_dir, ...
	  ana_dir, ...
	  svd_params.contrast_no, ...
	  svd_params.eigen_name, ...
	  'default', ...
	  '', ...
	  svd_params.temporal_filter, ...
	  svd_params.normalize, ...
	  svd_params.project_residual);
  fclose(fid);
  MM('SVD', argfile);
  spm_unlink(argfile);
  
end
