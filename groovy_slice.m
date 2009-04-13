function groovy_slice(glob_ps, sub_ps)
% slice timing metabatch file
  
for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  s_filter = [glob_ps.prefixes.slice this_sub.raw_filter];

  % Reference slice 
  ref_slice = this_sub.slice.ref_slice;
  
  first_flag = 1;
  for ss = 1:length(this_sub.sesses) % and session 
    dirn = fullfile(glob_ps.fdata_root, ...
		    this_sub.dir, this_sub.sesses(ss).dir);
    P = spm_get('files', dirn, s_filter);
    
    % get information from first file, for first session
    if first_flag % first session for each subject
      first_img = deblank(P(1,:));
      V = spm_vol(first_img);
      
      % Sets slice time information
      % value 1 is time to acquire one slice
      % value 2 is time between beginning of last slice
      % and beginning of first slice of next volume
      st = this_sub.slice.time;
      sl_times = [st st + (this_sub.TR-st*V.dim(3))];
      first_flag = 0;
    end

    % do slice timing correction
    spm_slice_timing(P,...
		     this_sub.slice.acq_order,...
		     ref_slice,sl_times);

  end
end
  
