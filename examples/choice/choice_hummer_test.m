%  Test batch file for later transfer to the Hummer

% This path is the same on Hummer and Hoover
parameter_root = '/import/hoover/imagers/choice';

% You might want to change this to something like /var/tmp/choice for
% final testing, and for use on the Hummer itself
fdata_root = '/home/imagers/choice';

subject_sdir = '05AM';
batch_sdir = 'groove_1';
ana_sdir = 'spm2_test_ana';

addpath(fullfile(parameter_root, batch_sdir));

choice_one_subject;
