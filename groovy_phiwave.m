function groovy_phiwave(glob_ps, sub_ps)
% phiwave metabatch file
  
% store path
pwd_orig = pwd;

% get phiwave stuff from global parameters
params = glob_ps.phiwave;

for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);

  % get, goto SPM results directory
  sub_dir = fullfile(glob_ps.fdata_root,this_sub.dir);
  ana_sdir = fullfile(sub_dir, glob_ps.stats.ana_sdir);
  cd(ana_sdir);
  
  % Make phido design object
  pD = phido('SPM.mat');

  % Turn off masking?
  % (if params.no_mask is present, and non-zero)
  if groovy_struct('getifthere', params, 'no_mask')
    nScan = n_time_points(pD);
    xM = -ones(nScan,1)/0;  % fancy eh?  Not ones(nScan, 1) * -Inf note
    pD = masking_struct(pD, xM);
  end
  
  % Remove 's' prefix from image names in design, so we can work on the
  % unsmoothed images (which assumes they exist)
  pD = prefix_images(pD, 'remove', 's');
  
  % Estimate the design, doing wavelet transform on the way
  pE = estimate(pD, [], params.estimate);

  % Merge old contrasts
  xCon = get_contrasts(pD);
  pE = add_contrasts(pE, xCon(2:end));
  xCon = get_contrasts(pE);

  % Write new denoised images with default file name 
  % (standard deviation image would be prepended with 'std_')
  Ic = params.contrasts.numbers;
  for c = 1:length(Ic)
    c_i = Ic(c);
    c_name = mars_utils('str2fname', ...
			[xCon(c_i).name params.contrasts.suffix]);
    get_wdimg(pE, c_i, params.denoise, c_name);
  end
end
cd(pwd_orig);










