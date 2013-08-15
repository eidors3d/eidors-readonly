% Create 3D model of a Rectangular resistor
% $Id$

scale = .56; %meters
ll=0:scale:5;
ww=( 0:scale:2 ).^1.2;
hh=( 0:scale:3 ).^0.8;
conduc= .13;  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2; % Ohms

mdl= mk_grid_model([],ww,hh,ll);
mdl= rmfield(mdl,'coarse2fine');
mdl.boundary= find_boundary(mdl.elems);
mdl.gnd_node = 1;
elec_nodes= [1:length(ww)*length(hh)];
elec(1).nodes= elec_nodes;                  elec(1).z_contact= z_contact;
elec(2).nodes= num_nodes(mdl)-elec_nodes+1; elec(2).z_contact= z_contact;
stim.stim_pattern= [-1;1]*current;
stim.meas_pattern= [-1,1];
mdl.stimulation= stim;
mdl.electrode= elec;
show_fem(mdl);

% analytical solution
Len = max(ll) - min(ll);
Wid = max(ww) - min(ww);
Hig = max(hh) - min(hh);
Block_R =  Len / ( Wid * Hig ) / conduc;
Contact_R = z_contact;
R = Block_R + 2*Contact_R;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);


mdl.solve = @fwd_solve_1st_order;
mdl.system_mat = @system_mat_1st_order;

img= mk_image( mdl, conduc );

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
