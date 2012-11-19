% Forward Model of a cylindrical object
% Create 3D model of a tunnel $Id: tunnelsim01.m 2356 2010-11-08 08:41:42Z aadler $
n_elec = 32;
shape_str = ['solid incyl  = cylinder (0,0,0; 1,0,0; 1) -maxh=1.0; \n', ...
    'solid farcyl = cylinder (0,0,0; 1,0,0; 5) -maxh=5.0; \n' ...
    'solid pl1    =  plane(-5,0,0;-1,0,0);\n' ...
    'solid pl2    =  plane(5,0,0; 1,0,0);\n' ...
    'solid mainobj= pl1 and pl2 and farcyl and not incyl;\n'];
th= linspace(0,2*pi,n_elec+1)'; th(end)=[];
cth= cos(th); sth=sin(th); zth= zeros(size(th));
elec_pos = [zth, cth, sth, zth cth, sth];
elec_shape= 0.01;
elec_obj = 'incyl';
fmdl = ng_mk_gen_models(shape_str, elec_pos, elec_shape, elec_obj);


% Construct the Wenner stimulation pattern
fmdl.stimulation= stim_pattern_geophys( n_elec, 'Wenner', {'spacings', 1:6,'circumferential_meas',1} );

% Compute the geometrical factor for the apparent resistivity estimation
img1= mk_image(fmdl,1);
vh1= fwd_solve(img1);
normalisation= 1./vh1.meas;
I= speye(length(normalisation));
I(1:size(I,1)+1:size(I,1)*size(I,1))= normalisation;

% Construct a model with a homogeneous conductivity of 0.1 Sm and insert a conductive sphere
img = mk_image(fmdl,0+ mk_c2f_circ_mapping(fmdl,[0;1.5;0;0.3])*100);
img.elem_data(img.elem_data==0)= 0.1;
show_fem(img,3);
print_convert stimPatt_PseudoSec07_1.png

% Solve the forward problem
dd  = fwd_solve(img);

% Show the pseudo-section of the apparent resistivity
show_pseudosection( fmdl, I*dd.meas, {'CircularOutside','yz'})
print_convert stimPatt_PseudoSec07_2.png
