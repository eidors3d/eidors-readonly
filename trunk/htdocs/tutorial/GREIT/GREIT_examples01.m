% Simulate obj $Id$

fmdl = ng_mk_cyl_models([2,1,0.1],[16,1],[0.05]); 
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);

% Homogeneous object
imgh= mk_image( fmdl, 1 );
% Fwd Simulation
vh = fwd_solve( imgh);

% Target Object
select_fcn = inline('(x-0.2).^2 + (y-0.5).^2 + (z-1).^2<0.1^2','x','y','z');
memb_frac = elem_select( imgh.fwd_model, select_fcn);
imgt= mk_image( fmdl, 1 + memb_frac );
% Fwd Simulation
vt = fwd_solve( imgt);

% Add SNR 2.0 noise
vn = add_noise(4, vt, vh);

show_fem(imgt);
print_convert GREIT_examples01a.png
