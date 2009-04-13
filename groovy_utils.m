function varargout = groovy_utils(action, varargin)
% Perform various utility tasks for groovy batch system
%
% Input
% action    - string defining action to perform (see below)
% 
% Possible actions are:
% SPM = groovy_utils('make_ss_ttest', P)
% Makes a single sample ttest RFX SPM design
% Where:
% P        - string list of image names, or array of vol structs
% SPM      - returned single sample ttest SPM structure
% 
% [C, names] = groovy_utils('movement_regressors, rparam_file, mr_types)
% returns functions of movement parameters, and regressor names 
% for single SPM session, for later inclusion into SPM design 
% Where:
% rparam_file    - name of text file containing movement parameters for
%                  this SPM session
% mr_types       - defines what function of movement parameters you
%                  want.  Cell array.  Can contain none or more of
%                  'moves'  - just simple movement parameters 
%                  'mm1'    - movement parameters offset by one
%                  (i.e. row number r contains movement parameters
%                  relating scan r-1 to first scan in series)
%                  'moves_2' - movement parameters squared (.^2)
%                  'mm1_2'   - mm1 parameters squared
%                  Default value for mr_types is {'moves'}
%                  If you pass {} as mr_types, you will get empty 'C',
%                  and 'names' output matrices.
% 
% Output 
% C              - movement parameter regressors
% names          - cell array of names for regressors columns (one value
%                  per column in C matrix
% 
% $Id: groovy_utils.m,v 1.2 2005/12/30 11:52:35 matthewbrett Exp $ 
  
if nargin < 1
  error('Need action argument to groovy_utils')
end

switch lower(action)
 case 'make_ss_ttest'
  varargout{1} = sf_make_ss_ttest(varargin{:});
 case 'movement_regressors'
  [A B] = sf_movement_regressors(varargin{:});
  varargout = {A, B};
 otherwise
  error('Oh no, not that')
end

return
 
% Subfunctions which do the actual work

function SPM = sf_make_ss_ttest(P)
% makes default single sample ttest design from image list
%
% Input
% P       - list of images (strings or cell struct)
% 
% Returns
% SPM     - SPM structure
%
% $Id: groovy_utils.m,v 1.2 2005/12/30 11:52:35 matthewbrett Exp $ 
 
if nargin < 1
  P = spm_get(Inf, 'IMAGE', 'Select images for ttest');
end
if iscell(P)
  P = char(P);
end
if ischar(P)
  V = spm_vol(P); 
  P = cellstr(P);
elseif isstruct(P)
  V = P;
  P = {V(:).fname};
else
  error('Was expecting string or vol struct as input');
end

nscan = length(V);
SPM.xY.P = P;
SPM.xY.VY = V;
SPM.nscan = nscan;
SPM.xX = struct(...
    'X', ones(nscan, 1), ...
    'iH', 1, ...
    'iC', [], ...
    'iB', [], ...
    'iG', [], ...
    'name', {{'mean'}}, ...
    'I', [1:nscan; ones(3, nscan)]', ...
    'sF', {{'obs'  ''  ''  ''}} ...
    );
SPM.xC = [];
SPM.xGX = struct(...
    'iGXcalc',  1, ...
    'sGXcalc',  'omit', ...
    'rg',  [], ...
    'iGMsca',  9, ...
    'sGMsca',  '<no grand Mean scaling>', ...
    'GM',  0, ...
    'gSF',  ones(nscan, 1), ...
    'iGC',  12, ...
    'sGC',  '(redundant: not doing AnCova)', ...
    'gc',  [], ...
    'iGloNorm',  9, ...
    'sGloNorm',  '<no global normalisation>'...
    );
SPM.xVi = struct(...
    'iid', 1, ...
    'V', speye(nscan));
SPM.xM = struct(...
    'T', -Inf, ...
    'TH', ones(nscan, 1) * -Inf, ...
    'I', 1, ...
    'VM', [], ...
    'xs', struct(...
	'Analysis_threshold', 'None (-Inf)', ...
	'Implicit_masking', 'Yes: NaNs treated as missing', ...
	'Explicit_masking', 'No')...
    );
P3 = sprintf('leaving %d degrees of freedom from %d images', ...
	     nscan - 1, nscan);
SPM.xsDes = struct(...
    'Design',  'One sample t-test', ...
    'Global_calculation',  'omit', ...
    'Grand_mean_scaling',  '<no grand Mean scaling>', ...
    'Global_normalisation',  '<no global normalisation>', ...
    'Parameters', {{'1 condition, +0 covariate, +0 block, +0 nuisance', ...
		   '1 total, having 1 degrees of freedom', ...
		   P3}});
SPM.SPMid = 'SPM2: spm_spm_ui (v2.49)';

return

function [C, names] = sf_movement_regressors(rparam_file, mr_types)
% returns functions of movement parameters, and regressor names 

if nargin < 1
  rparam_file = spm_get([0 1], 'rp*.txt', ...
			'Select movement parameter file');
end
if nargin < 2
  mr_types = {'moves'};
end

% Names for movement parameters
move_names = {'x trans', 'y trans', 'z trans', ...
	       'x rot', 'y rot', 'z rot'};
cmn = char(move_names);

moves = spm_load(rparam_file);
if isempty(moves)
  error(['Cannot get movement parameters from: ' rparam_file]);
end
if size(moves, 2) ~= 6
  error(['Expecting a Nx6 movement parameter matrix from file' ...
	rparam_file]);
end

% mean centre moves
moves = moves - ones(size(moves, 1), 1) * mean(moves);
moves_m1 = [zeros(1, 6); moves(1:end-1,:)];

C = []; 
names = {}; 
    
if ismember('moves', mr_types)
  C = [C moves];
  names = [names move_names];
end
if ismember('mm1', mr_types)
  C = [C moves_m1];
  names = [names sf_name_cat(move_names, ' minus 1')];
end
if ismember('moves_2', mr_types)
  C = [C moves.^2];
  names = [names sf_name_cat(move_names, ' .^2')];
end
if ismember('mm1_2', mr_types)
  C = [C moves_m1.^2];
  names = [names sf_name_cat(move_names, ' minus 1 .^2')];
end
return

function names = sf_name_cat(n, suffix)
names = cellstr([char(n) repmat(suffix, 6, 1)])';
return

