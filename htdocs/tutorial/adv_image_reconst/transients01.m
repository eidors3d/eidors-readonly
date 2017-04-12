%Create 2D model
ball1 = 'solid ball = cylinder(0.1,0.4,0;0.1,0.4,1;0.3) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;';
box1  = 'solid box = orthobrick(-0.6,-0.6,0; 0.1,-0.1,0.05) -maxh=0.1;';
fmdl= ng_mk_cyl_models(0,[8],[0.2,0,0.05],{'ball','box',[ball1,box1]}); 
fmdl.stimulation = stim_meas_list([1,5,2,3],8,.01,1);
inclusion = vertcat( fmdl.mat_idx{2:3} );
%Conductivity and permittivity parameters
img = mk_image( fmdl, 1);

img.elem_data(inclusion) = .01;
img.fwd_solve.get_all_meas=1; %Get all measurements

subplot(221); show_fem(img);
print_convert transients01a.png
