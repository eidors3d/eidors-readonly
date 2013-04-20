function data= np_fwd_solve( fwd_model, img)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% Input:
%    fwd_model = forward model
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

warning('EIDORS:deprecated','NP_FWD_SOLVE is deprecated as of 07-Jun-2012. Use FWD_SOLVE_1ST_ORDER instead.');

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

Vfwd = forward_solver(s_mat.E, p.I, tol, s_mat.perm);

Velec=Vfwd( p.n_node+(1:p.n_elec),:);
voltH = zeros( p.n_meas, 1 );
idx=0;
for i=1:p.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   voltH( idx+(1:n_meas) ) = meas_pat*Velec(:,i);
   idx= idx+ n_meas;
end

% create a data structure to return
data.meas= voltH;
data.time= NaN; % unknown
data.name= 'solved by np_fwd_solve';
try; if img.fwd_solve.get_all_meas == 1
   data.volt = Vfwd(1:p.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = Vfwd;               % all, including CEM nodes
end; end

function do_unit_test
   test_3d_resistor

function test_3d_resistor
   ll=5; % length
   ww=1; % width
   hh=1; % height
   conduc= .13;  % conductivity in Ohm-meters
   current= 4;  % Amps
   z_contact= 9e-1;
   scale = .46;
   nn=0;
   for z=0:ll; for x=0:ww; for y=0:hh
      nn=nn+1;
      mdl.nodes(nn,:) = [x,y,z];
   end; end; end

   mdl= mk_grid_model([],0:ww,0:hh,0:ll);
   mdl.nodes= mdl.nodes*scale;
   mdl= rmfield(mdl,'coarse2fine');
   mdl.boundary= find_boundary(mdl.elems);
   mdl.gnd_node = 1;
   mdl.normalize_measurements = 0;
   elec_nodes= [1:(ww+1)*(hh+1)];
   elec(1).nodes= elec_nodes;      elec(1).z_contact= z_contact;
   elec(2).nodes= nn-elec_nodes+1; elec(2).z_contact= z_contact;
   stim.stim_pattern= [-1;1]*current;
   stim.meas_pattern= [-1,1];
   mdl.stimulation= stim;
   mdl.electrode= elec;
%  show_fem(mdl);

   % analytical solution
   R = ll / ww / hh / scale/ conduc + 2*z_contact/scale^2;

   V= current*R;
%  fprintf('Solver %s: %f\n', 'analytic', V);

   mdl.solve = @np_fwd_solve;
   mdl.system_mat = @np_calc_system_mat;
%  mdl.misc.perm_sym= '{n}'; %%% Check it works without
   img= mk_image(mdl, conduc);

   fsol= fwd_solve(img);
%  fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('np_fwd_solve vs analytic', fsol.meas, V, 1e-11);



   mdl.solve = @fwd_solve_1st_order;
   mdl.system_mat = @system_mat_1st_order;

   img= mk_image(mdl, conduc);

   fsol= fwd_solve(img);
%  fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('new solver 2d vs analytic', fsol.meas, V, 1e-11);
