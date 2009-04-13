% Single subject batch script
%
% Assumes the following parameters (here with example contents)
% parameter_root = '/some/directory/on/my/machine';
% node_root = '/another/directory/on/my/machine';
% subject_sdir = 'subject_1';
% batch_sdir = 'my_batch_dir';
% ana_sdir = 'my_ana_sdir';

% Need to add groovy_batch to path if not already present
% addpath /usr/local/spm/groovy_batch

% Get the global and subject parameters
glob_ps.parameter_root = parameter_root;
glob_ps.fdata_root = fdata_root;
[glob_ps sub_ps] = choice_top_groove(subject_sdir, glob_ps);
glob_ps.ana_sdir = ana_sdir;

groovy_reorient(glob_ps, sub_ps);
groovy_realign(glob_ps, sub_ps);
groovy_coreg(glob_ps, sub_ps);
groovy_norm_calc(glob_ps, sub_ps);
groovy_norm_write(glob_ps, sub_ps);
groovy_smooth(glob_ps, sub_ps);
groovy_subject_model(glob_ps, sub_ps);
groovy_contrasts(glob_ps, sub_ps)
