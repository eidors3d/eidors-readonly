% Construct the Schlumberger stimulation pattern
spacing= [1 1 1 2 3 4 6 8 8 11 12 14 17]; 		% Spacing between electrodes (usually called the "a" parameter)
multiples= [1 2 3 2 5/3 6/4 7/6 1 10/8 1 13/12 15/14 1];% Multiples (usually called the "n" parameter)
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Schlumberger', {'spacings', spacing,'multiples',multiples} );

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[100;0;-50;50])*100);
img.elem_data(img.elem_data==0)= 0.1;

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( fmdl, dd.meas);
print_convert stimPatt_PseudoSec03.png

