function groovy_subject_model(glob_ps, sub_ps, estimateyn)
% metabatch file to run SPM2 FMRI models for single subjects
if nargin < 3
  estimateyn = 1;
end
defs = glob_ps.defaults.stats.fmri;  
  
% store path
pwd_orig = pwd;
   
% Specify design stuff common to all subjects
%===========================================================================
% global normalization: OPTOINS:'Scaling'|'None'
%---------------------------------------------------------------------------
gSPM.xGX.iGXcalc       = 'None';

% low frequency confound: high-pass cutoff (secs) [Inf = no filtering]
%---------------------------------------------------------------------------
gSPM.xX.K(1).HParam    = glob_ps.stats.high_pass;

% intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w'
%-----------------------------------------------------------------------
%
%SPM.xVi.form       = 'AR(1) + w';
gSPM.xVi.form       = glob_ps.stats.ar_form;

% basis functions and timing parameters
%---------------------------------------------------------------------------
% OPTIONS:'hrf'
%         'hrf (with time derivative)'
%         'hrf (with time and dispersion derivatives)'
%         'Fourier set'
%         'Fourier set (Hanning)'
%         'Gamma functions'
%         'Finite Impulse Response'
%---------------------------------------------------------------------------
% Fill in the field below with the corresponding string above
gSPM.xBF.name       = glob_ps.stats.event_bf.name;

% length in seconds - not used for hrf 
gSPM.xBF.length     = glob_ps.stats.event_bf.length; 

% order of basis set - not used for hrf
gSPM.xBF.order      = glob_ps.stats.event_bf.order;  

% The next two fields usually don't need changing. 
% number of time bins per scan
gSPM.xBF.T          = 16;                % number of time bins per scan
gSPM.xBF.T0         = 1;                 % first time bin (see slice timing) - middle of TA

% Selfish explanatory - OPTIONS: 'scans'|'secs' for onsets
gSPM.xBF.UNITS      = glob_ps.stats.units;

% Order of convolution OPTIONS: 1|2  
% 1 means no Volterra.  This is the right answer.
gSPM.xBF.Volterra   = 1;                 

for sb = 1:length(sub_ps)
  this_sub = sub_ps(sb);
  
  % Get template SPM structure
  SPM = gSPM;
  
  % Note that the TR must be the same for all runs in a model
  SPM.xY.RT          =  this_sub.TR;    % seconds
  
  % specify filter for filenames
  Filter             = ['^' glob_ps.prefixes.stats this_sub.raw_filter];
  
  % get, make, goto SPM results directory
  sub_dir = fullfile(glob_ps.fdata_root, this_sub.dir);
  ana_sdir = fullfile(sub_dir, glob_ps.stats.ana_sdir);
  if ~(exist(ana_sdir))
    mkdir(sub_dir,glob_ps.stats.ana_sdir);
  end
  cd(ana_sdir);
  
  PP=[];
  
  for ss = 1:length(this_sub.sesses)
    % Information for this session
    this_ss = this_sub.sesses(ss);
    
    % directory containing scans
    fildir = fullfile(sub_dir, this_ss.dir);
    
    % file selection
    P = spm_select('FPList', fildir, Filter);
    n_scans = size(P,1);
    SPM.nscan(ss) = n_scans;
    PP = strvcat(PP, P);
    
    % Condition stuff - onsets, durations, types.
    for cno = 1:length(this_ss.cond_names)
      ons = this_ss.ons{cno};
      dur = this_ss.dur{cno};
      parameters = this_ss.parameters{cno};
      if isempty(parameters)
	P = struct('name', 'none');
      else
	P = struct('name', 'other',...
		   'P', parameters(:), ...
		   'h', 1); % h is order, 1 being just linear
    end

      early_ons = ons < 0;
      if any(early_ons)
	warning('Some onsets are less than 0');
	disp(ons(early_ons))
      end
      late_ons = ons >= n_scans;
      if any(late_ons), 
	warning('Some onsets are after the end of scanning');
	disp(ons(late_ons))
      end
      ons = ons(~(early_ons | late_ons));
      if isempty(ons)
	error(...
	    sprintf('Session %d, %d scans, cond %d; no onsets are left',...
		ss, n_scans,cno));
      end
      dur = dur(~(early_ons | late_ons));
      SPM.Sess(ss).U(cno) = struct(...
	  'ons', ons, ...
	  'dur', dur,...
	  'name',{this_ss.cond_names(cno)},...
	  'P', P); % Parametric modulation
    end
    
    % Movement stuff
    movefil = spm_select('FPList', fildir, '^rp.*\.txt$');
    if isempty(movefil)
      error(['Cannot find movement parameters with ' movefil]);
    end
    [m_c m_names] = groovy_utils('movement_regressors', ...
				 movefil, ...
				 glob_ps.stats.movement_params);
    
    if ~isempty(m_c)
      % Fix any goofy realignment params
      m_c = m_c(1:n_scans, :);
    end
    
    % other covariates and names go in first
    covs = [this_ss.covs m_c];
    cov_names = [this_ss.cov_names m_names];
    
    % Put covariates into model
    SPM.Sess(ss).C.C    = covs;     % [n x c double] covariates
    SPM.Sess(ss).C.name = cov_names; % [1 x c cell]   names
    
  end % session loop
  
  % set files
  SPM.xY.P           = PP;
  
  % Configure design matrix
  SPM = spm_fmri_spm_ui(SPM);

  % Set explicit mask if required
  if groovy_struct('isthere', glob_ps.stats, 'explicit_mask')
    nimgs = size(SPM.xY.VY);
    imgv = ones(nimgs);
    SPM.xM.T =  imgv * -Inf;
    SPM.xM.TH = imgv * -Inf;
    SPM.xM.I  = 0;
    
    msk = glob_ps.stats.explicit_mask;
    SPM.xM.VM = spm_vol(msk);
  end
  
  if estimateyn
    % Estimate parameters
    spm_unlink(fullfile('.', 'mask.img')); % avoid overwrite dialog
    SPM = spm_spm(SPM);
  else
    save SPM SPM
  end
end

% back to initial directory
cd(pwd_orig);
