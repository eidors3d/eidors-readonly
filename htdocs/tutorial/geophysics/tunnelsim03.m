% Create 3D model of a tunnel $Id$ 

% Reconstruct entire image (slow)
% Simulation protocol. TODO: we need a geophysical stim protocol
imdl = mk_common_model('d2c2',N_elec);
imdl.fwd_model = fmdl;
imdl.jacobian_bkgnd.value = cond_mdl;

imgr = inv_solve( imdl, vs_h, vs_i );

imgr.calc_colours.npoints= 128;
slices = [0.0,inf,inf,1,1;
          0.5,inf,inf,2,1; 
          1.0,inf,inf,3,1]; 
subplot(211); show_slices(imgr,slices);
print_convert tunnelsim03a.png
