% Forward Model including a borehole

borehole_diameter= 0.1;
n_elec= 32;
e0 = -linspace(0.5,5,n_elec)';
elec_pos = [0*e0,0*e0+borehole_diameter,e0,0*e0,1+0*e0,0*e0];
elec_shape= [0.04 0.02 0.1];

incyl = sprintf('cylinder (0,0,0; 0,0,-6;%f)',borehole_diameter); % Borehole
farcyl= sprintf('cylinder (0,0,0; 0,0,-15;10)'); % Medium surrounding the borehole
shape_str = ['solid incyl  = ',incyl,' -maxh=0.1 ; \n', ...
             'solid farcyl = ',farcyl,' -maxh=10 ; \n', ...
             'solid pli1   =  plane(0,0,0;0,0,1);\n' ...
             'solid pli2   =  plane(0,0,-6; 0,0,-1);\n' ...
             'solid plu1   =  plane(0,0,0;0,0,1);\n' ...
             'solid plu2   =  plane(0,0,-15; 0,0,-1);\n' ...
             'solid innerobj= pli1 and pli2 and incyl;\n', ...
             'solid mainobj= plu1 and plu2 and farcyl and not innerobj;\n'];
clear elec_obj
for i=1:size(e0,1) ;  elec_obj{i} = 'incyl';end
[fmdl,mat_idx] = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);


% Construct the Wenner stimulation pattern
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:25} );

% Use apparent_resistivity
fmdl.jacobian = @jacobian_apparent_resistivity;
fmdl.solve    = @fwd_solve_apparent_resistivity;


% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[0;0.7;-3;0.3])*100);
img.elem_data(img.elem_data==0)= 0.1;
show_fem(img);
print_convert stimPatt_PseudoSec05_1.png '-density 75'
 
% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
fmdl.show_pseudosection.orientation = 'vertical';
show_pseudosection( fmdl, dd.meas);
print_convert stimPatt_PseudoSec05_2.png '-density 75'
