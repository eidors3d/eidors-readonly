function test_2d_resistor(opt)
% Create 2D model of a cylindrical resistor
% $Id$

if nargin>0 && strcmp(opt,'UNIT_TEST'); do_unit_test; return; end

resistor_test;

function img = this_mdl(conduc, z_contact, current);
   nn= 12;     % number of nodes
   ww=3;       % width = 4
   scale = .35;
   mdl=mk_grid_model([],3+scale*(1:ww), scale*(1:nn/ww));
   mdl= rmfield(mdl,'coarse2fine'); % don't calc this.

   mdl.gnd_node = 1;
   elec_nodes= [1:ww];
   elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
   elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
   stim.stim_pattern= [-1;1]*current;
   stim.meas_pattern= [-1,1];
   mdl.stimulation= stim;
   mdl.electrode= elec;
   n_el = size(mdl.elems,1);
   img= eidors_obj('image','2D rectangle', ...
         'elem_data', ones(n_el,1) * conduc );
   img.fwd_model = mdl;
   img.fwd_model.normalize_measurements = 0;

function V = analytic_soln(img);
   % analytical solution
   nodes = img.fwd_model.nodes;
   wid_len= max(nodes) - min(nodes);
   conduc =  mean(img.elem_data);
   Block_R = wid_len(2) / wid_len(1) / conduc;
   % Contact R reflects z_contact / width. There is no need to scale
   %  by the scale, since this is already reflected in the size of the
   %  FEM as created by the grid. This is different to the test_3d_resistor,
   %  where the FEM is created first, and then scaled, so that the ww
   %  and hh need to be scaled by the scale parameter.
   z_contact = sum([img.fwd_model.electrode(:).z_contact]);
   current = max(img.fwd_model.stimulation(1).stim_pattern(:));
   Contact_R = z_contact/wid_len(1);
   R = Block_R + Contact_R;

   V= current*R;
   fprintf('Solver %s: %f\n', 'analytic', V);
   fprintf('Solver %s: %f\n', 'analytic (no z_contact)', V - 2*Contact_R*current);

function V = first_order_solver( img );
   % AA_SOLVER
   img.fwd_model.solve = @fwd_solve_1st_order;
   img.fwd_model.system_mat = @system_mat_1st_order;
   fsol= fwd_solve(img);
   fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   V = fsol.meas;

function vals = resistor_test;

   current= 4;  % Amps
   conduc=  .4; % conductivity in Ohm-meters
   z_contact= 1e-1;
   img = this_mdl( conduc, z_contact, current);

   show_fem(img.fwd_model);

   vals.analytic = analytic_soln(img);

% AA_SOLVER
   vals.aa_solver = first_order_solver(img);

warn_state = warning('query','EIDORS:Deprecated');
warning('off','EIDORS:Deprecated');


% NP_SOLVER
img.fwd_model.solve = @np_fwd_solve;
img.fwd_model.system_mat = @np_calc_system_mat;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
vals.np_solver = fsol.meas;

% NOW CALCULATE THE ANALYTICAL JACOBIAN
img.fwd_model.solve = @np_fwd_solve;
img.fwd_model.jacobian = @np_calc_jacobian;

img.elem_data= ones(num_elems(img),1) * conduc ;
Jnp= calc_jacobian(img);

img.fwd_model.jacobian = @jacobian_perturb;
Jp1= calc_jacobian(img);

img.elem_data= ones(num_elems(img),1) * conduc ;
Jp2= zeros(size(Jnp));
delta= 1e-8;
for i=1:num_elems(img)
   fsol_h= fwd_solve(img);
   img.elem_data(i) = img.elem_data(i) + delta;
   fsol_i= fwd_solve(img);
   Jp2(:,i) = (fsol_i.meas - fsol_h.meas)/delta; 
end
   
mdl.solve = @fwd_solve_1st_order;
mdl.jacobian = @jacobian_adjoint;
Jaa= calc_jacobian(img);


fprintf('Jacobians: Cols by Jaa, Jnp, Jp1, Jp2:\n')
JJ = [Jaa;Jnp;Jp1;Jp2];
disp(JJ(:,1:6))




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
mdl.solve = @fwd_solve_1st_order;
mdl.system_mat = @system_mat_1st_order;
img.fwd_model = mdl;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
% NP_SOLVER
mdl.solve = @np_fwd_solve;
mdl.system_mat = @np_calc_system_mat;
fsol= fwd_solve(img);
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
fprintf('Analytic is not expected to be same in last case\n');

warning(warn_state.state,'EIDORS:Deprecated');


function do_unit_test
  vals = resistor_test;

  unit_test_cmp('Analytic vs AA', vals.analytic, vals.aa_solver, 1e-10);
  unit_test_cmp('Analytic vs NP', vals.analytic, vals.np_solver, 1e-10);

