function groovy_fieldmap(glob_ps, sub_ps)
% fieldmap metabatch file
  
for sb = 1:length(sub_ps) % for each subject
  this_sub = sub_ps(sb);
  epi_img = this_sub.fieldmap.target;
  
  params = groovy_struct('merge', glob_ps.fieldmap, this_sub.fieldmap);
  
  prev_imgs = [];
  % Run the fieldmap stuff - can be different for different sessions
  for ss = 1:length(this_sub.sesses)
    this_sess  = this_sub.sesses(ss);
    params = groovy_struct('merge', params, this_sess.fieldmap);
    
    % Get images to work on
    these_imgs = params.imgs;
    
    % Check if images are different from previous session, don't rerun if
    % not
    do_f = any(size(these_imgs) ~= size(prev_imgs));
    if ~do_f
      do_f = ~all(all(these_imgs == prev_imgs)); 
    end
    if ~do_f, continue, end
    
    sf_process_fieldmap(params, these_imgs, epi_img);
    prev_imgs = these_imgs;
    
  end
end
return

function sf_process_fieldmap(pm_defs, fm_imgs, epi_img)
% function to process fieldmaps
  
if nargin < 3
  error('Need defs, field map and epi image');
end

IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Define parameters for fieldmap creation
IP.et{1} = pm_defs.SHORT_ECHO_TIME;
IP.et{2} = pm_defs.LONG_ECHO_TIME;

% Set parameters for unwrapping
IP.uflags.iformat = pm_defs.INPUT_DATA_FORMAT;
IP.uflags.method = pm_defs.UNWRAPPING_METHOD;
IP.uflags.fwhm = pm_defs.FWHM;
IP.uflags.pad = pm_defs.PAD;
IP.uflags.ws = pm_defs.WS;
IP.uflags.etd = pm_defs.LONG_ECHO_TIME - pm_defs.SHORT_ECHO_TIME;     

% Set parameters for unwarping 
IP.ajm = pm_defs.DO_JACOBIAN_MODULATION;
IP.blipdir = pm_defs.K_SPACE_TRAVERSAL_BLIP_DIR;
IP.tert = pm_defs.TOTAL_EPI_READOUT_TIME;
IP.epifm = pm_defs.EPI_BASED_FIELDMAPS;

% Clear any old handles etc
IP.fm = [];
IP.vdm = [];
IP.jim = [];
IP.pP = [];
IP.epiP = [];
IP.uepiP = [];
IP.vdmP = [];
ID = cell(4,1);

%----------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
%----------------------------------------------------------------------

if ischar(fm_imgs), fm_imgs = spm_vol(fm_imgs); end
n_fms = length(fm_imgs);
switch n_fms
 case 4  % real, imaginary pairs
  for i = 1:n_fms, IP.P{i} = fm_imgs(i); end

  %----------------------------------------------------------------------
  % Create field map (in Hz) - this routine calls the unwrapping
  %----------------------------------------------------------------------

  IP.fm = FieldMap('CreateFieldMap',IP);
  
  %----------------------------------------------------------------------
  % Write out field map
  % Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
  %----------------------------------------------------------------------
  
  FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Smoothed phase map');

 case 2  % precalculated Hz map
  IP.pP = fm_imgs(1);
  IP.fmagP = fm_imgs(2);
  IP.fm.fpm = spm_read_vols(IP.pP);
  IP.fm.jac = pm_diff(IP.fm.fpm,2);
  if isfield(IP,'P') & ~isempty(IP.P{1})
    IP.P = cell(1,4);
  end
  otherwise 
   error('Funny number of input fieldmap images')
end

%----------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

[IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);

%----------------------------------------------------------------------
% Select an EPI to unwarp
%----------------------------------------------------------------------

if ischar(epi_img), epi_img = spm_vol(epi_img); end
IP.epiP = epi_img;

%----------------------------------------------------------------------
% Match voxel displacement map to image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

if isfield(pm_defs, 'match_vdm')
  if pm_defs.match_vdm
    IP.vdmP = FieldMap('MatchVDM',IP);
  end
end

%----------------------------------------------------------------------
% Unwarp EPI
%----------------------------------------------------------------------

IP.uepiP = FieldMap('UnwarpEPI',IP);

%----------------------------------------------------------------------
% Write unwarped EPI 
% Outputs -> uNAME-OF-EPI.img
%----------------------------------------------------------------------
IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dim(4),'Unwarped image');

return









