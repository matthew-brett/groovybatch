% Metametabatch file to run all FIAC analysis
% (see http://phiwave.sourceforge.net/fiac)
% You need:
% The unwarp toolbox in the path ahead of SPM 
% Also on path:
% FieldMap toolbox (with FieldMap.m replaced by FieldMap.m at
%    http://phiwave.sourceforge.net/fiac/FieldMap.m)
% MarsBaR (marsbar.sourceforge.net) (for Phiwave) - v0.38.2
% Phiwave (phiwave.sourceforge.net) v3.3
% Groovy_batch (http://phiwave.sourceforge.net/groovy_batch) v0.1

% At first we start using all the subjects' data (block and event)
[g s] = fiac_top_groove;
groovy_segment(g, s);
groovy_slice(g, s);
groovy_realign(g, s);
groovy_fieldmap(g, s);
groovy_unwarp(g, s);
groovy_norm_calc(g, s);
groovy_norm_write(g, s);

% Now just the block experiments
[g s] = fiac_top_groove('','',struct('experiment_types', 'block'));

% Different smoothing for single and rfx analysis 
ss_smooth = 5;
ss_sdir = 'ana_5mm';
rfx_smooth = 10;
rfx_sdir = 'ana_10mm';

% All phiwave analyses in rfx sdirs
phi_sdir = rfx_sdir;

% Do the example single subject, with 5mm smoothing
g.stats.FWHM = ss_smooth;
g.stats.ana_sdir = ss_sdir;
groovy_smooth(g, s(1));
groovy_subject_model(g, s(1));
groovy_contrasts(g, s(1));

% And the whole lot at 10mm smoothing for random effects analysis
g.stats.FWHM = rfx_smooth;
g.stats.ana_sdir = rfx_sdir;
groovy_smooth(g, s);
groovy_subject_model(g, s);
groovy_contrasts(g, s);
groovy_randeff(g, s);

% Run phiwave model for all subjects (on unsmoothed images)
groovy_phiwave(g, s);

% (on my machine this gives a mkdir error sometimes; resolved by just
% rerunning script from here)
groovy_phiwave_randeff(g, s);

% Finally you might want to run the figures scripts
fiac_apply_masking(g, s, phi_sdir);
fiac_write_spm_std(g, s, ss_sdir, rfx_sdir);
fiac_figures(g, s, ss_sdir, rfx_sdir);
