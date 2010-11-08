% Create 3D model of a tunnel $Id$ 

% Reconstruct entire image (slow)
% Simulation protocol. TODO: we need a geophysical stim protocol
imdl = mk_common_model('d2c2',N_elec);
imdl.fwd_model = fmdl;
imdl.jacobian_bkgnd.value = cond_mdl;

imgr = inv_solve( imdl, vs_h, vs_i );

show_fem(imgr); ylim(2*[-1,1]); zlim(2*[-1,1]);

view(90,0); print_convert tunnelsim03a.png
