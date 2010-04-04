% RPI tank model $Id$

% simple inverse model -> replace fields to match this model
imdl = mk_common_model('b2c2',32);
imdl.fwd_model.normalize_measurements = 0;

imdl.fwd_model.electrode = imdl.fwd_model.electrode([8:-1:1, 32:-1:9]);

imdl.fwd_model = rmfield(imdl.fwd_model,'meas_select');
imdl.fwd_model.stimulation = stim;
imdl.hyperparameter.value = 1;


load Rensselaer_EIT_Phantom;
vh =  real(ACT2006_homog); vi = real(ACT2000_phant);

img = inv_solve(imdl, vh, vi);
img.calc_colours.cb_shrink_move = [0.5,0.8,0];
show_fem(img,[1,1]); axis off; axis image

print_convert('rpi_data02a.png','-density 60');
