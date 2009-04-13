function P = fiac_write_spm_std(g, s, ss_spm_sdir, rfx_spm_sdir)
% makes error (std) images for selected SPM contrasts 

fiac_root = g.fdata_root;
t_filt = 'spmT_*.img';
prefix = 'de_t_';

con_names = {s(1).contrasts.name};

P = spm_get('files', ...
	    fullfile(fiac_root, s(1).dir, ss_spm_sdir), ...
	    t_filt);

for dn = 1:length(con_names)
  sdir = sprintf('rfx_%s_%s', rfx_spm_sdir, ...
	      mars_utils('str2fname', con_names{dn}));
  PP = spm_get('files', fullfile(fiac_root, sdir), t_filt);
  P = strvcat(P, PP);
end

for i = 1:size(P, 1)
  t_img = spm_vol(deblank(P(i,:)));
  [pn fn ext] = fileparts(t_img.fname);
  fn = strrep(fn, 'spmT_', 'con_');
  con_img = spm_vol(fullfile(pn, [fn ext]));
  
  cimg = spm_read_vols(con_img);
  img = spm_read_vols(t_img);
  img = 1 / (img ./ cimg);
  oimg = t_img;
  [pn fn ext] = fileparts(oimg.fname);
  oimg.fname = fullfile(pn, [prefix fn ext]);
  spm_write_vol(oimg, img);
end

