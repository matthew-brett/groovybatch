function SPM = set_globals_ar(SPM, global_calc, ar_form)
% Set global thresholding and AR options for SPM
% FORMAT SPM = set_globals_ar(SPM, ar_form)  
% This is copy / paste from spm_fmri_spm_ui
  
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model fillup',0);
spm_help('!ContextHelp',mfilename)
cVi = ar_form;
K = SPM.xX.K;
nscan = [];
nsess = numel(K);
for sno = 1:numel(K)
  nscan(sno) = length(K(sno).row);
end

% create Vi struct
if ~ischar(cVi)	% AR coeficient[s] specified
  SPM.xVi.Vi = spm_Ce(nscan,cVi(1:3));
  cVi        = ['AR( ' sprintf('%0.1f ',cVi) ')'];
else
  switch lower(cVi)
   case 'none'		%  xVi.V is i.i.d
    SPM.xVi.V  = speye(sum(nscan));
    cVi        = 'i.i.d';
   otherwise		% otherwise assume AR(0.2) in xVi.Vi
    SPM.xVi.Vi = spm_Ce(nscan,0.2);
    cVi        = 'AR(0.2)';
    end
end
SPM.xVi.form = cVi;


%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% get file identifiers
%=======================================================================

%-Map files
%-----------------------------------------------------------------------
fprintf('%-40s: ','Mapping files')                          	     %-#
VY    = spm_vol(SPM.xY.P);
fprintf('%30s\n','...done')                                 	     %-#


%-check internal consistency of images
%-----------------------------------------------------------------------
spm_check_orientations(VY);

%-place in xY
%-----------------------------------------------------------------------
SPM.xY.VY = VY;

%-Compute Global variate
%=======================================================================
GM    = 100;
q     = length(VY);
g     = zeros(q,1);
fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
for i = 1:q
  fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,q)) %-#
  g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')               %-#

% scale if specified (otherwise session specific grand mean scaling)
%-----------------------------------------------------------------------
gSF   = GM./g;
if strcmp(lower(global_calc),'none')
  for i = 1:nsess
    row = SPM.xX.K(i).row
    gSF(row) = GM./mean(g(row));
  end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for i = 1:q
	SPM.xY.VY(i).pinfo(1:2,:) = SPM.xY.VY(i).pinfo(1:2,:)*gSF(i);
end

%-place global variates in global structure
%-----------------------------------------------------------------------
SPM.xGX = struct(...
    'iGXcalc', global_calc,...
    'rg', g, ...
    'GM', GM,...
    'gSF', gSF);

%-Masking structure automatically set to 80% of mean
%=======================================================================
try
	TH    = g.*gSF*defaults.mask.thresh;
catch
	TH    = g.*gSF*0.8;
end
SPM.xM        = struct(	'T',	ones(q,1),...
			'TH',	TH,...
			'I',	0,...
			'VM',	{[]},...
			'xs',	struct('Masking','analysis threshold'));
