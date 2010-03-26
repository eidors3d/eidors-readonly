% Netgen simulation $Id$
subplot(221)

% Get a basic inverse model - replace the fwd model part
imdl = mk_common_model('a2c2',16); % the parameters aren't important because we replace it
imdl.fwd_model = img.fwd_model;

imgr = inv_solve(imdl, vh, vi);
show_fem(imgr);
axis equal; print -dpng -r100 netgen_sims06a.png

% Change the hyperparameter
imdl.hyperparameter.value = .003;
imgr = inv_solve(imdl, vh, vi);
show_fem(imgr);
axis equal; print -dpng -r100 netgen_sims06b.png
