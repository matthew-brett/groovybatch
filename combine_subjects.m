function nSPM = combine_subjects(varargin)
% Makes event model combined across sessions
% FORMAT nSPM = combine_sessions(SPM1, SPM2, ...)
%
% Inputs
% SPM1,SPM2 etc - spm FMRI design structures for subjects
%
% Output
% nSPM    - modified fMRI design combining across subjects
%
% Covariance estimation form and global calcalation method
% taken from first design
%
% Typical usage
% >> subj1 = load('subj1/SPM');
% >> subj2 = load('subj2/SPM');
% >> subj1 = load('subj3/SPM');
% >> mkdir combined_subs
% >> cd combined_subs
% >> SPM = combine_subjects(subj1.SPM, subj2.SPM, subj3.SPM);
% >> save SPM SPM
% >> SPM = spm_spm(SPM); % estimate design

if nargin < 2
  error('Need at least two SPMs to combine')
end
spms = varargin;
col_offset = 0;
row_offset = 0;

SPM1 = spms{1};
nxY = struct('RT', SPM1.xY.RT, 'P', [], 'VY', []);
nxX = struct('X', [],...
	     'iB', [],...
	     'iG', [], ...
	     'name', [],...
	     'K', []);
nscan = [];
Sess = [];
for spmno = 1:numel(spms)
  SPM = spms{spmno};
  if ~isstruct(SPM)
    error(sprintf('SPM number %d is not a structure', spmno))
  end
  % Relevant stuff from this SPM
  xX = SPM.xX;
  xY = SPM.xY;
  % Make design matrix
  X = xX.X;
  [m n] = size(X);
  nxX.X = blkdiag(nxX.X, X);
  % Get, fix filter row numbers
  K = xX.K;
  n_sess = numel(K);
  for nk = 1:numel(K)
    row = K(nk).row;
    new_row = row + row_offset;
    K(nk).row = new_row;
    Sess(end+1).row = new_row;
    nscan(end+1) = length(row);
  end
  nxX.K = [nxX.K K];
  % Set block effect and covariate columns
  nxX.iB = [nxX.iB xX.iB+col_offset];
  nxX.iG = [nxX.iG xX.iG+col_offset];
  % Names
  nxX.name = [nxX.name xX.name];
  % Add data
  if nxY.RT ~= xY.RT
    error('RTs do not match across sessions');
  end
  nxY.P = strvcat(nxY.P, xY.P);
  % Move to next subject row set
  row_offset = row_offset + m;
  col_offset = col_offset + n;
end
nSPM = struct(...
    'SPMid', SPM1.SPMid,...
    'swd', pwd, ...
    'xY', nxY,...
    'xX', nxX,...
    'Sess', Sess,...
    'nscan', nscan);
nSPM = set_globals_ar(nSPM, 'session specific', SPM1.xVi.form);
