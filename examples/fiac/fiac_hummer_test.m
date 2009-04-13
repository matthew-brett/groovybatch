%  Test batch file for later transfer to the Hummer
parameter_root = '/import/hoover/imagers/FIAC';
fdata_root = '/home/imagers/FIAC';
subject_sdir = 'fiac3';
batch_sdir = 'too_gru';
ana_sdir = 'spm2_test';

addpath(fullfile(parameter_root, batch_sdir));

fiac_one_subject;
