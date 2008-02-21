% Create 2D model of a cylindrical resistor
% $Id: test_2d_resistor.m,v 1.3 2008-02-21 19:06:50 aadler Exp $

nn= 12;     % number of nodes
ww=2;       % width = 4
conduc= .1;  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2;
mdl= eidors_obj('fwd_model','2D rectangle');
mdl.nodes = [floor( (0:nn-1)/ww );rem(0:nn-1,ww)]';
mdl.elems = delaunayn(mdl.nodes);
mdl.boundary= find_boundary(mdl.elems);
mdl.gnd_node = 1;
elec_nodes= [1:ww];
elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
stim.stim_pattern= [-1;1]*current;
stim.meas_pattern= [-1,1];
mdl.stimulation= stim;
mdl.electrode= elec;
show_fem(mdl);
img= eidors_obj('image','2D rectangle', ...
      'elem_data', ones(size(mdl.elems,1),1) * conduc );

% AA_SOLVER
mdl.solve = @aa_fwd_solve;
mdl.system_mat = @aa_calc_system_mat;
img.fwd_model = mdl;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% AA_SOLVER
mdl.solve = @np_fwd_solve;
mdl.system_mat = @np_calc_system_mat;
img.fwd_model = mdl;
fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% analytical solution
wid_len= max(mdl.nodes) - min(mdl.nodes);
R = wid_len(1) / wid_len(2) / conduc + 2*z_contact;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);
