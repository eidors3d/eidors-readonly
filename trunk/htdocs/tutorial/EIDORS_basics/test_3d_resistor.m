% Create 3D model of a Rectangular resistor
% $Id$

ll=5*1; % length (# elements)
ww=1*2; % width (# elements)
hh=1*3; % height (# elements)
conduc= .13;  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2; % Ohms
scale = .46; %

mdl= mk_grid_model([],0:ww,0:hh,0:ll);
mdl= rmfield(mdl,'coarse2fine');
mdl.nodes= mdl.nodes*scale;
mdl.boundary= find_boundary(mdl.elems);
mdl.gnd_node = 1;
elec_nodes= [1:(ww+1)*(hh+1)];
elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
stim.stim_pattern= [-1;1]*current;
stim.meas_pattern= [-1,1];
mdl.stimulation= stim;
mdl.electrode= elec;
show_fem(mdl);

% analytical solution
Block_R =  ll / ww / hh / scale/ conduc;
Contact_R = z_contact;
R = Block_R + 2*Contact_R;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);


mdl.solve = @fwd_solve_1st_order;
mdl.system_mat = @system_mat_1st_order;

img= mk_image( mdl, conduc );

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
