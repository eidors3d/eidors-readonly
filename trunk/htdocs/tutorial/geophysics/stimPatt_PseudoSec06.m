% Forward Model of a cylindrical object
n_elec= 32;
ring_vert_pos = [0.5]; 
fmdl= ng_mk_cyl_models([1,0.3,0.05],[n_elec,ring_vert_pos],[0.02,0.05,0.02]);


% Construct the Wenner stimulation pattern
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:7,'circumferential_meas',1} );

% Compute the geometrical factor for the apparent resistivity estimation
img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[0;0.05;0.5;0.1])*100);
img.elem_data(img.elem_data==0)= 0.1;
show_fem(img,3);
print_convert stimPatt_PseudoSec06_1.png
 
% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( fmdl, I*dd.meas, 'CircularInside')
print_convert stimPatt_PseudoSec06_2.png

