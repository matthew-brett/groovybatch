#!/usr/bin/python
#
# Script to fetch archives, unpack to temporary directory, set
# up matlab batch run, and copy results back to main file space

import os, sys, shutil

server_root = '/import/hoover/imagers/choice/archives'
parameter_root = '/import/hoover/imagers/choice'
fdata_root = '/var/tmp/choice'
one_sub_mfile = 'choice_one_subject'

# Name of the analysis subdirectory for SPM
ana_sdir = 'spm2_humm_ana'

# Name of batch directory in parameter_root
batch_sdir = 'groove_1'

# We need the list of all subjects to work out which subject to use here.
subjects = (
    '03FR',
    '04AL',
    '05AM',
    '06AW',
    '07CS',
    '09NM',
    '10PL',
    '11KS',
    '12KK',
    '13AM',
    '14AD',
    '15BD',
    '16NK',
#   '17HI',   - massive artefacts - subject excluded
    '18RL'
        )

# Task ID from environment variable
task_id = os.environ.get('SGE_TASK_ID')
if not task_id:
    task_id = 1

# Get the subject number from the SGE task ID
task_no = int(task_id)
subject = subjects[task_no-1]

# make the fdata_root directory if it doesn't exist
# then go there.
subj_dir = "%s/%s" % (fdata_root, subject)
try:
    os.makedirs(subj_dir)
except OSError:
    pass
os.chdir(fdata_root)

# Download and unpack the bz2'ed data
tar_file =  "%s.bz2" % (subject,)  # tar file for this subject
shutil.copy("%s/%s" % (server_root, tar_file), '.')
os.system("tar jxvf %s" % tar_file)
os.unlink(tar_file)

# cd to the analysis (NOT BATCH)  directory
ana_dir_full = "%s/%s" % (subj_dir, ana_sdir)
try:
    os.mkdir(ana_dir_full)
except OSError:
    pass
os.chdir(ana_dir_full)

# make matlab startup batch file, and run batch job
start_file = "startup.m"
f = open(start_file, 'wt')
f.write("""%% Python generated matlab startup file
parameter_root = '%s';
fdata_root = '%s';
subno = %d
subject_sdir = '%s';
ana_sdir = '%s';
batch_sdir = '%s'; 

addpath(fullfile(parameter_root, batch_sdir));
%s;

exit
""" % (parameter_root, fdata_root, task_no, subject,
       ana_sdir, batch_sdir, one_sub_mfile))
f.close()
os.system('matlab -nojvm')
os.remove(start_file)

# Pack model directory into tar archive, and go
os.chdir(fdata_root)
ana_tar = "%s_%s.tar.gz" % (subject, ana_sdir) 
os.system("tar zcvf %s %s/%s" % (ana_tar, subject, ana_sdir))
shutil.move(ana_tar, server_root)

# Delete analysis directory to save space
shutil.rmtree(subject, ignore_errors=True)
