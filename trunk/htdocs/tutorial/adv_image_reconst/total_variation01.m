% Create mesh with blocky objects $Id$

% Simulation (forward model) 1024 elements
imdl= mk_common_model('d2c0',16); 
fmdl= imdl.fwd_model;

% Identify block in centre
ctrs= interp_mesh(fmdl);
xe= ctrs(:,1); ye= ctrs(:,2);
re= sqrt(xe.^2+ye.^2);
block=(ye>0 & re<.69); % for 1024

sim_img= mk_image(fmdl,1);
v_homg= fwd_solve(sim_img);

sim_img.elem_data(block) = 1.1;
v_simu= fwd_solve(sim_img);

clf;
calc_colours('greylev',-.1); % white background
show_fem(sim_img)
print -r50 -dpng total_variation01a.png;
