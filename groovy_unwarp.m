function groovy_unwarp(glob_ps, sub_ps)
% realignment / unwarp metabatch file
  
% Get realignment, unwarp defaults
defs = glob_ps.defaults.realign;
unwarp = glob_ps.defaults.unwarp;

% Flags to pass to routine to calculate realignment parameters
% (spm_realign)
reaFlags = struct(...
    'Quality', defs.estimate.quality,...  % estimation quality
    'fwhm', 5,...                         % smooth before calculation;
    'rtm', 0,...                          % whether to realign to mean 
    'PW',''...                            % weighting images for each subject
    );

% Flags to pass to routine to create resliced images
% (spm_reslice)
resFlags = struct(...
    'interp', defs.write.interp,...       % trilinear interpolation
    'wrap', defs.write.wrap,...           % wrapping info (ignore...)
    'mask', defs.write.mask,...           % masking (see spm_reslice)
    'which', glob_ps.realign.which, ...   % images to write
    'mean', glob_ps.realign.write_mean);  % whether to write mean image

% Unwarp options parsing
uwe_flags = struct('order',    unwarp.estimate.basfcn,...
		   'sfP',      [],...
		   'regorder',    unwarp.estimate.regorder,...
		   'lambda',      unwarp.estimate.regwgt,...
		   'jm',          unwarp.estimate.jm,...
		   'fot',         glob_ps.fieldmap.foe,...
		   'sot',         glob_ps.fieldmap.soe,...
		   'fwhm',        unwarp.estimate.fwhm,...
		   'rem',         unwarp.estimate.rem,...
		   'noi',         unwarp.estimate.noi,...
		   'exp_round',   unwarp.estimate.expround);

uwr_flags = struct('interp',      defs.write.interp,...
		   'wrap',        defs.write.wrap,...
		   'mask',        defs.write.mask,...
		   'which',       resFlags.which,...
		   'mean',        resFlags.mean);

if unwarp.estimate.jm == 1
  uwr_flags.udc = 2;
else
  uwr_flags.udc = 1;
end

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  n_sesses = length(this_sub.sesses);
  pp  = cell(1, n_sesses);
  ppm = pp;
  r_filter = [glob_ps.prefixes.realign this_sub.raw_filter];
  for ss = 1:n_sesses % and session 
    this_sess = this_sub.sesses(ss);
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sess.dir);
    p = spm_get('files', dirn, r_filter);
    pp{ss}  = p;
    if isfield(this_sess, 'vdm_img')
      if ~isempty(this_sess.vdm_img)
	ppm{ss} = this_sess.vdm_img;
      end
    end
  end

  if glob_ps.fieldmap.do_realign
    spm_realign(pp, reaFlags);
  end

  clear ads
  tmpP = spm_vol(pp{1}(1,:));
  uwe_flags.M = tmpP.mat;
  for ss=1:n_sesses
    uwe_flags.sfP = ppm{ss};
    ds = spm_uw_estimate(pp{ss}, uwe_flags);
    ads(ss) = ds;
    [path,name,ext,ver] = fileparts(pp{ss}(1,:));
    pefile = fullfile(path,[name '_uw.mat']);
    save(pefile,'ds');
  end		
  
  % And do reslice
  spm_uw_apply(ads, uwr_flags);
  
end











