% Construct the Wenner stimulation pattern
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:32} );

% Use apparent_resistivity
fmdl.jacobian = @jacobian_apparent_resistivity;
fmdl.solve    = @fwd_solve_apparent_resistivity;

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[100;0;-50;50])*100);
img.elem_data(img.elem_data==0)= 0.1;

show_fem(img,3);
print_convert stimPatt_PseudoSec01.png

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( fmdl, dd.meas);
 print_convert stimPatt_PseudoSec02.png

