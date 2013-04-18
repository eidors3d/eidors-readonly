% Forward Model of a cylindrical object
n_elec= 32;
ring_vert_pos = [0.5]; 
fmdl= ng_mk_cyl_models([1,0.3,0.05],[n_elec,ring_vert_pos],[0.02,0.05,0.02]);


% Construct the Wenner stimulation pattern
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:7,'circumferential_meas',1} );

% Use apparent_resistivity
fmdl.jacobian = @jacobian_apparent_resistivity;
fmdl.solve    = @fwd_solve_apparent_resistivity;

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0.1+ mk_c2f_circ_mapping(fmdl,[0;0.10;0.5;0.1])*100);
show_fem(img);
print_convert stimPatt_PseudoSec06_1.png '-density 75'
 
% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
fmdl.show_pseudosection.orientation = 'CircularInside';
show_pseudosection( fmdl, dd.meas);
print_convert stimPatt_PseudoSec06_2.png '-density 75'
