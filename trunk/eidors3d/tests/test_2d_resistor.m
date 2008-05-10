% Create 2D model of a cylindrical resistor
% $Id: test_2d_resistor.m,v 1.6 2008-05-10 17:42:38 aadler Exp $

nn= 12;     % number of nodes
ww=2;       % width = 4
if ~exist('conduc');conduc=   1;end  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2;
scale = .35;
mdl=mk_grid_model([],3+scale*(1:ww), scale*(1:nn/ww));
mdl.gnd_node = 1;
elec_nodes= [1:ww];
elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
stim.stim_pattern= [-1;1]*current;
stim.meas_pattern= [-1,1];
mdl.stimulation= stim;
mdl.electrode= elec;
show_fem(mdl);
n_el = size(mdl.elems,1);
img= eidors_obj('image','2D rectangle', ...
      'elem_data', ones(n_el,1) * conduc );

% AA_SOLVER
mdl.solve = @aa_fwd_solve;
mdl.system_mat = @aa_calc_system_mat;
img.fwd_model = mdl;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% NP_SOLVER
mdl.solve = @np_fwd_solve;
mdl.system_mat = @np_calc_system_mat;
img.fwd_model = mdl;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% analytical solution
wid_len= max(mdl.nodes) - min(mdl.nodes);
R = wid_len(1) / wid_len(2) / conduc + 2*z_contact/scale;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);

% NOW CALCULATE THE ANALYTICAL JACOBIAN
mdl.solve = @np_fwd_solve;
mdl.jacobian = @np_calc_jacobian;
mg.elem_data= ones(n_el,1) * conduc ;
Jnp= calc_jacobian(mdl,img);

mdl.jacobian = @perturb_jacobian;
Jp1= calc_jacobian(mdl,img);

img.elem_data= ones(n_el,1) * conduc ;
Jp2= zeros(size(Jnp));
delta= 1e-8;
for i=1:n_el
   fsol_h= fwd_solve(img);
   img.elem_data(i) = img.elem_data(i) + delta;
   fsol_i= fwd_solve(img);
   Jp2(:,i) = (fsol_i.meas - fsol_h.meas)/delta; 
end
   
mdl.solve = @aa_fwd_solve;
mdl.jacobian = @aa_calc_jacobian;
Jaa= calc_jacobian(mdl,img);

[Jaa;Jnp;Jp1;Jp2]




% TEST RESISTOR THAT IS NOT RECTANGULAR
s1=.6;s2=3*s1;
m2= eidors_obj('fwd_model','not rect');
m2.nodes= [s1/2,-1;-s1/2,-1;s1/2,0;-s1/2,0;s2/2,1;-s2/2,1];
m2.elems= [1,2,3;2,3,4;3,4,5;4,5,6];
m2.gnd_node= 3;
m2.electrode(1).nodes= [1,2];
m2.electrode(2).nodes= [3,4];
m2.electrode(3).nodes= [5,6];
[m2.electrode(:).z_contact]= deal(.0001);
m2.stimulation(1).stim_pattern= [-1;0;1];
m2.stimulation(1).meas_pattern= [1,-1,0;0,-1,1];
% AA_SOLVER
m2.solve = @aa_fwd_solve;
m2.system_mat = @aa_calc_system_mat;
img= eidors_obj('image','','fwd_model',m2,'elem_data',[1,1,1,1]);
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
