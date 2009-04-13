% Single subject batch script
%
% Assumes the following parameters (here with example contents)
% parameter_root = '/some/directory/on/my/machine';
% node_root = '/another/directory/on/my/machine';
% subject_sdir = 'subject_1';
% batch_sdir = 'my_batch_dir';
% ana_sdir = 'my_ana_sdir';

addpath /usr/local/spm/too_groovy
addpath /usr/local/spm/unwarp
addpath /usr/local/spm/spm2/toolbox/FieldMap
addpath /usr/local/spm/marsbar
addpath /usr/local/spm/phiwave

% Get the global and subject parameters
g.parameter_root = parameter_root;
g.fdata_root = fdata_root;
[g s] = fiac_top_groove(subject_sdir, g);
g.stats.ana_sdir = ana_sdir;

groovy_segment(g, s);
groovy_slice(g, s);
groovy_realign(g, s);
groovy_fieldmap(g, s);
groovy_unwarp(g, s);
groovy_norm_calc(g, s);
groovy_norm_write(g, s);
groovy_smooth(g, s);

% Now just the block experiment
just_block = struct('experiment_types', 'block');
[g s] = fiac_top_groove(subject_sdir,g, just_block);
g.stats.ana_sdir = ana_sdir;
groovy_subject_model(g, s);
groovy_contrasts(g, s);
groovy_phiwave(g, s);

