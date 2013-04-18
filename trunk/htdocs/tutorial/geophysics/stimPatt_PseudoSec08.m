% Forward Model of a cylindrical object
% Create 3D model of a tunnel $Id$
n_elec = 32;
shape_str = ['solid incyl  = cylinder (0,0,0; 1,0,0; 1) -maxh=1.0; \n', ...
    'solid pl1    =  plane(-.5,0,0;-1,0,0);\n' ...
    'solid pl2    =  plane( .5,0,0; 1,0,0);\n' ...
    'solid mainobj= pl1 and pl2 and incyl;\n'];
th= linspace(0,2*pi,n_elec+1)'; th(end)=[];
cth= cos(th); sth=sin(th); zth= zeros(size(th));
elec_pos = [zth, cth, sth, zth cth, sth];
elec_shape= 0.01;
elec_obj = 'incyl';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);


% Construct the Wenner stimulation pattern
 fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings',1:6,'circumferential_meas',1} );
%fmdl.stimulation= mk_stim_patterns(n_elec, 1, [0,3],[0,3],{},1);

% Use apparent_resistivity
fmdl.jacobian = @jacobian_apparent_resistivity;
fmdl.solve    = @fwd_solve_apparent_resistivity;

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0.1+ mk_c2f_circ_mapping(fmdl,[0;0.8;0;0.1])*100);
show_fem(img);
print_convert stimPatt_PseudoSec08_1.png '-density 75'

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
fmdl.show_pseudosection.orientation = {'CircularInside','yz'};
show_pseudosection( fmdl, dd.meas);

print_convert stimPatt_PseudoSec08_2.png '-density 75'
