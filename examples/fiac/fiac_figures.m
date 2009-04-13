function fiac_figures(g, s, ss_spm_sdir, ss_phi_sdir)
% Script to make figures for FIAC phiwave analysis

global MF
MF.fiac_root = g.fdata_root;
MF.ss_spm_dir = fullfile(MF.fiac_root, s(1).dir, ss_spm_sdir);
MF.ss_phi_dir = fullfile(MF.fiac_root, s(1).dir, ss_phi_sdir);
MF.rfx_dir_pattern = ['rfx_' ss_phi_sdir '_%s'];

% Contrast names 
MF.con_names = {'Effects of interest', s(1).contrasts.name};

% Image types
MF.types = {'effect', 'error', 't'};

% Colormaps for effect, error, t
MF.types_cmaps = {'flow.lut', hot(256), 'flow.lut'};

% Antomicals
MF.ss_anat.vol = spm_vol(fullfile(MF.fiac_root, 'fiac1', 'anat1', ...
			       'wfiac1_anat1.img'));

MF.ss_anat.range = [0 20000];
MF.grp_anat.vol = spm_vol(fullfile(spm('dir'), 'canonical', ...
				'avg152T1.mnc'));
MF.grp_anat.range = [0 1];

% Single subject thing
range_fudge=0.6;

so = slover;
so.figure = figure;
so.img(1).cmap = gray;
so.printstr = 'print -dpng -noui';
so.img(2).cmap = 'flow.lut';
so.cbar = 2;

level_labels = {'s1', 'rfx'};
st_labels = {'phi', 'spm'};

% For figures 1 and 2
% 1) FSt contrast, S1, SPM vs Phiwave, Effect, error, t, axial slice z=4
% 2) FSt contrast, RFX, SPM vs Phiwave, Effect, error, t, axial slice z=4
so.slices = 4;
cname = 'FSt';
fcname = sf_str2fname(cname);
for level = [1:2] % and for figures 1:2
  [phi_st spm_st anat_st] = sf_get_vols(cname, ...
					level-1, ...
					[0 1 1], ...
					'centile', ...
					0.1);
  sts = [phi_st spm_st];
  for t = 1:length(MF.types)
    for s = 1:length(sts)
      so = sf_fill_fields(so, anat_st, sts(s).(MF.types{t}));
      so = paint(so);
      g_name = sprintf('%s_%s_%s_%s.png', ...
		       fcname, ...
		       level_labels{level}, ...
		       MF.types{t}, ...
		       st_labels{s});
      print_fig(so, g_name);
    end
  end
end

% 3) a)  DSt-SSt, RFX, SPM vs Phiwave, Effect, axial slice z=4
level = 2; % rfx
t = 1; % effect only
so.slices = 4;
cname = 'Dst-SSt';
fcname = sf_str2fname(cname);
[phi_st spm_st anat_st] = sf_get_vols(cname, ...
				      level-1, ...
				      [0 1 1], ...
				      'centile', ...
				      0.1);
sts = [phi_st spm_st];
for s = 1:length(sts)
  so = sf_fill_fields(so, anat_st, sts(s).(MF.types{t}));
  so = paint(so);
  g_name = sprintf('%s_%s_%s_%s.png', ...
		   fcname, ...
		   level_labels{level}, ...
		   MF.types{t}, ...
		   st_labels{s});
  print_fig(so, g_name);
end

% 3) b)  DSp-SSp, RFX, SPM vs Phiwave, Effect, axial slice z=-4
so.slices = -4;
cname = 'DSp-SSp';
fcname = sf_str2fname(cname);

[phi_st spm_st anat_st] = sf_get_vols(cname, ...
				      level-1, ...
				      [0 1 1], ...
				      'centile', ...
				      0.1);
sts = [phi_st spm_st];
for s = 1:length(sts)
  so = sf_fill_fields(so, anat_st, sts(s).(MF.types{t}));
  so = paint(so);
  g_name = sprintf('%s_%s_%s_%s.png', ...
		   fcname, ...
		   level_labels{level}, ...
		   MF.types{t}, ...
		   st_labels{s});
  print_fig(so, g_name);
end

% 4) Interaction - RFX, SPM vs Phiwave, Effect, for:
% Coronal Y = -60, Sagittal X = 0, Axial Z = -14.
cname = 'Interaction';
fcname = sf_str2fname(cname);
[phi_st spm_st anat_st] = sf_get_vols(cname, ...
				      level-1, ...
				      [0 1 1], ...
				      'centile', ...
				      0.1);
sts = [phi_st spm_st];
slices = [-60 0 -10];
orientations = {'coronal', 'sagittal', 'axial'};
for o = 1:length(orientations)
  for s = 1:length(sts)
    so = sf_fill_fields(so, anat_st, sts(s).(MF.types{t}));
    so.transform = orientations{o};
    so.slicedef = [];
    so.slices = slices(o);
    so = paint(so);
    g_name = sprintf('%s_%s_%s_%s.png', ...
		     fcname, ...
		     orientations{o}, ...
		     MF.types{t}, ...
		     st_labels{s});
    print_fig(so, g_name);
  end
end

return

% Use hackey global variables to get the right stuff
function [phi_st, spm_st, anat_st] = sf_get_vols(cname, rfx_f, match_range,calc_type, arg)					
global MF
if nargin < 3
  match_range = zeros(1, length(MF.types));
end

if ~rfx_f
  con_num = find(ismember(MF.con_names, cname));
  if isempty(con_num), error('Cannot find contrast'); end
  spm_st.effect.vol.fname = fullfile(MF.ss_spm_dir, sprintf('con_%.4d.img', con_num));
  phi_st.effect.vol.fname = fullfile(MF.ss_phi_dir, sprintf('masked_%s_stein.img', cname));
  spm_st.error.vol.fname = fullfile(MF.ss_spm_dir, sprintf('de_t_spmT_%.4d.img', con_num));
  phi_st.error.vol.fname = fullfile(MF.ss_phi_dir, sprintf('masked_std_%s_stein.img', cname));
  spm_st.t.vol.fname = fullfile(MF.ss_spm_dir, sprintf('spmT_%.4d.img', con_num));
  phi_st.t.vol.fname = fullfile(MF.ss_phi_dir, sprintf('masked_t_%s_stein.img', cname));
  anat_st = MF.ss_anat;
else
  fcname = sf_str2fname(cname);
  con_num = 2;
  cname = 'mean';
  spm_dir = fullfile(MF.fiac_root, sprintf(MF.rfx_dir_pattern, fcname));
  phi_dir = fullfile(MF.fiac_root, sprintf('rfx_wv_%s_stein', fcname));
  spm_st.effect.vol.fname = fullfile(spm_dir, sprintf('con_%.4d.img', con_num));
  phi_st.effect.vol.fname = fullfile(phi_dir, sprintf('masked_%s_stein.img', cname));
  spm_st.error.vol.fname = fullfile(spm_dir, sprintf('de_t_spmT_%.4d.img', con_num));
  phi_st.error.vol.fname = fullfile(phi_dir, sprintf('masked_std_%s_stein.img', cname));
  spm_st.t.vol.fname = fullfile(spm_dir, sprintf('spmT_%.4d.img', con_num));
  phi_st.t.vol.fname = fullfile(phi_dir, sprintf('masked_t_%s_stein.img', cname));  
  anat_st = MF.grp_anat;
end
for t = 1:length(MF.types)
  et = MF.types{t};
  spm_st.(et).vol = spm_vol(spm_st.(et).vol.fname);
  phi_st.(et).vol = spm_vol(phi_st.(et).vol.fname);
  if match_range(t)
    rng = sf_get_range(calc_type, arg, spm_st.(et).vol, phi_st.(et).vol);
    spm_st.(et).range = rng;
    phi_st.(et).range = rng;
  else
    spm_st.(et).range = sf_get_range(calc_type, arg, spm_st.(et).vol);
    phi_st.(et).range = sf_get_range(calc_type, arg, phi_st.(et).vol);
  end
  if strcmp(et, 'error')
    spm_st.(et).range(1) = 0;
    phi_st.(et).range(1) = 0;
  end
  spm_st.(et).cmap = MF.types_cmaps{t};
  phi_st.(et).cmap = MF.types_cmaps{t};
end

return

function rng = sf_get_range(type, arg, varargin)
mxmx = -Inf;
for i = 1:length(varargin)
  switch type
   case 'mxmn'
    [mx mn] = slover('volmaxmin', varargin{i});
   case 'centile'
    [mx mn] = sf_cent(varargin{i}, arg);
   otherwise
    error('Funny haha type')
  end
  mxmx = max([mxmx abs([mx mn])]);
end
rng = [-mxmx mxmx];
return

function fname = sf_str2fname(str)
% forbidden chars in file name
badchars = unique([filesep '/\ :;.''"~*?<>|&']);

tmp = find(ismember(str, badchars));   
if ~isempty(tmp)
  str(tmp) = '_';
  dt = diff(tmp);
  if ~isempty(dt)
    str(tmp(dt==1))=[];
  end
end
fname = str;
return

function so = sf_fill_fields(so, varargin)
for v = 1:length(varargin)
  st = varargin{v};
  fns = fieldnames(st);
  for f = 1:length(fns)
    fn = fns{f};
    so.img(v).(fn) = st.(fn);
  end
end
return

function [hi, lo] = sf_cent(vol, centile)
img = spm_read_vols(vol);
img = img(isfinite(img) & img~=0);
img = sort(img(:));
L = length(img);
pts = round(L .* ([centile 100-centile]/100));
hi = img(pts(2));
lo = img(pts(1));
return
