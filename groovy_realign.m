function groovy_realign(glob_ps, sub_ps)
% realignment metabatch file
  
% Get realignment defaults
defs = glob_ps.defaults.realign;

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

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  r_filter = [glob_ps.prefixes.realign this_sub.raw_filter];
  clear imgs; 
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_get('files', dirn, r_filter);
    imgs(ss) = {P};
  end
  
  % Run the realignment
  spm_realign(imgs, reaFlags);
  
  % Run the reslicing
  spm_reslice(imgs, resFlags);
end











