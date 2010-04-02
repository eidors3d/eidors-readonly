% RPI tank model $Id$

% simple inverse model -> replace fields to match this model
imdls = mk_common_model('b2c2',32);
imdls.fwd_model.normalize_measurements = 0;
imdls.fwd_model = rmfield(imdls.fwd_model,'meas_select');
imdls.fwd_model.stimulation = stim;
imdls.hyperparameter.value = .08;

imdl = imdls; imdl.fwd_model = fmdl; 

load Rensselaer_EIT_Phantom;
vh =  real(ACT2006(2:end)); vi = real(ACT2000(2:end));

subplot(221);
img = inv_solve(imdl , vh, vi);
show_fem(img);

%print -dpng -r125 rpi_data01a.png
print -depsc2  jnk.eps;!LD_LIBRARY_PATH="" convert -density 125 jnk.eps rpi_data02a.png

subplot(221);
img = inv_solve(imdls, vh, vi);
show_fem(img);
