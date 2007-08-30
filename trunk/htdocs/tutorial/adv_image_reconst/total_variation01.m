% Create mesh with blocky objects $Id: total_variation01.m,v 1.3 2007-08-30 03:58:27 aadler Exp $

% Simulation (forward model) 1024 elements
imdl= mk_common_model('d2c0',16); 
fmdl= imdl.fwd_model;

fm_nelems= size(fmdl.elems,1);
% Identify block in centre
xe= mean(reshape(fmdl.nodes(fmdl.elems',1),3,[]));
ye= mean(reshape(fmdl.nodes(fmdl.elems',2),3,[]));
re= sqrt(xe.^2+ye.^2);
block=(ye>0 & re<.69); % for 1024

sim_img= eidors_obj('image','sim_img','fwd_model',fmdl);
sim_img.elem_data= ones(fm_nelems,1);
v_homg= fwd_solve(sim_img);

sim_img.elem_data(block) = 1.1;
v_simu= fwd_solve(sim_img);

clf;
calc_colours('greylev',-.1); % black background
show_fem(sim_img)
print -r50 -dpng total_variation01a.png;
