% simulate radial movement $Id$
use_3d_model = 1;

if use_3d_model;
   fmdl= ng_mk_cyl_models([30,15,1],[16,15],[0.5]);
else;
   imdl = mk_common_model('f2d3c',16); fmdl= imdl.fwd_model;
end

stim_pat = mk_stim_patterns(16,1,'{ad}','{ad}', {'no_meas_current'}, 1);

fmdl.stimulation = stim_pat;
