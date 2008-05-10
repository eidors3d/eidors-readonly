% Create 2D model of a cylindrical resistor
% $Id: test_2d_resistor.m,v 1.7 2008-05-10 18:37:51 aadler Exp $

nn= 12;     % number of nodes
ww=2;       % width = 4
if ~exist('conduc');conduc=  .4;end  % conductivity in Ohm-meters
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
R = wid_len(2) / wid_len(1) / conduc + 2*z_contact/scale;

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
nn= 12;     % number of nodes
ww=2;       % width = 4
if ~exist('conduc');conduc=  .4;end  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2;
scale = .35;
mdl=mk_grid_model([],3+scale*(1:ww), scale*(1:nn/ww));
xdelta= .05*mdl.nodes(1:ww:nn,2);
mdl.nodes(1:ww:nn,1) = mdl.nodes(1:ww:nn,1) + xdelta;
mdl.nodes(2:ww:nn,1) = mdl.nodes(2:ww:nn,1) - xdelta;
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
img= eidors_obj('image','2D rectangle', 'fwd_model',mdl, ...
      'elem_data', ones(n_el,1) * conduc );
% AA_SOLVER
mdl.solve = @aa_fwd_solve;
mdl.system_mat = @aa_calc_system_mat;
fsol= fwd_solve(mdl,img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
% NP_SOLVER
mdl.solve = @np_fwd_solve;
mdl.system_mat = @np_calc_system_mat;
fsol= fwd_solve(mdl,img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% Analytic
e1= mdl.nodes(mdl.electrode(1).nodes,:);
d1= sqrt( sum(( [1,-1]*e1 ).^2 ));
e2= mdl.nodes(mdl.electrode(2).nodes,:);
d2= sqrt( sum(( [1,-1]*e2 ).^2 ));
hig= max(mdl.nodes(:,2)) - min(mdl.nodes(:,2));
% R = int_h (1/d) dh 
% d= a+b*h; a= d1; b= (d2-d1)/hig
% R = int 1/(a+b*h)*dh = log(a+bh)/b
% R = [log(d2) - log(d1)]*hig/(d2-d1)
% doesn't account for arc in conductivity pattern
R = (log(d2)-log(d1))*hig/(d2-d1)/conduc + 2*z_contact/scale;
V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);
