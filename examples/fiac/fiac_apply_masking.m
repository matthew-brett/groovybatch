function P = fiac_apply_masking(g, s, ss_sdir);
% makes masked version of stein images for display

fiac_root = g.fdata_root;
dnc_str = g.phiwave.denoise.thcalc;
wv_prefix = g.phiwave.estimate.wtprefix;
ss_dnc_filt = sprintf('*_%s.img', dnc_str);
rfx_dnc_filt = sprintf('*mean_%s.img', dnc_str);
V_msk = g.stats.explicit_mask;
prefix = 'masked_';

con_names = {s(1).contrasts.name};

P = spm_get('files', ...
	    fullfile(fiac_root, s(1).dir, ss_sdir), ...
	    ss_dnc_filt);

for dn = 1:length(con_names)
  dir = sprintf('rfx_%s%s_%s', ...
		wv_prefix, ...
		mars_utils('str2fname', con_names{dn}),...
		dnc_str);
  PP = spm_get('files', fullfile(fiac_root, dir), rfx_dnc_filt);
  P = strvcat(P, PP);
end

mimg = spm_read_vols(V_msk);
imgs = spm_vol(P);
for i_no = 1:length(imgs)
  img = spm_read_vols(imgs(i_no));
  img = img .* mimg;
  oimg = imgs(i_no);
  [pn fn ext] = fileparts(oimg.fname);
  oimg.fname = fullfile(pn, [prefix fn ext]);
  spm_write_vol(oimg, img);
end

