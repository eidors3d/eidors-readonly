% Create 3D model of a Rectangular resistor
% $Id$

ll=5*1; % length
ww=1*2; % width
hh=1*3; % height
conduc= .13;  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-1;
scale = .46;
mdl= eidors_obj('fwd_model','3D rectangle');
nn=0;
for z=0:ll; for x=0:ww; for y=0:hh
   nn=nn+1;
   mdl.nodes(nn,:) = [x,y,z];
end; end; end

if 0
% Matlab's delaunayn triangularization is screwed up
   mdl.elems = delaunayn(mdl.nodes);
elseif 0
% Matlab's delaunay3 triangularization is screwed up
   mdl.elems = delaunay3(mdl.nodes(:,1), mdl.nodes(:,2), mdl.nodes(:,3));
elseif 0
   elem1= [ 1 2 3 5; 5 6 2 3; 5 6 7 3; ...
            4 2 3 8; 8 6 2 3; 8 6 7 3;];
   mdl.elems=[];
   for i=0:ll-1;
       mdl.elems= [mdl.elems; elem1+4*i];
   end
else
   mdl= mk_grid_model([],0:ww,0:hh,0:ll);
   mdl.nodes= mdl.nodes*scale;
   mdl= rmfield(mdl,'coarse2fine');
end
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
mdl = mdl_normalize(mdl,0);



mdl.solve = @fwd_solve_1st_order;
mdl.system_mat = @system_mat_1st_order;

n_el = size(mdl.elems,1);
img= eidors_obj('image','3D rectangle', ...
      'elem_data', ones(n_el,1) * conduc, ...
      'fwd_model', mdl); 

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% analytical solution
Block_R =  ll / ww / hh / scale/ conduc;
Contact_R = z_contact/(ww*hh)/scale^2;
% Contact R reflects z_contact / (width/scale)^2. Here we need to use
%  the scale, since this is not reflected in the size of the
%  FEM as created by the grid. This is different to the test_2d_resistor,
%  where the FEM is created scaled, so that the ww
%  don't need to be scaled by the scale parameter.
R = Block_R + 2*Contact_R;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);
fprintf('Solver %s: %f\n', 'analytic (no z_contact)', V - 2*Contact_R*current);

if 0 % OLD CODE
   mdl.solve = @np_fwd_solve;
   mdl.jacobian = @np_calc_jacobian;
else
   mdl.solve = @fwd_solve_1st_order;
   mdl.jacobian= @jacobian_adjoint;
end
mdl.misc.perm_sym= '{n}';
mdl.normalize_measurements = 0;
img= eidors_obj('image','3D rectangle', ...
      'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
      'fwd_model', mdl); 

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);


% NOW CALCULATE THE ANALYTICAL JACOBIAN
if 0 % OLD CODE
   mdl.solve = @np_fwd_solve;
   mdl.jacobian = @np_calc_jacobian;
else
   mdl.solve = @fwd_solve_1st_order;
   mdl.jacobian= @jacobian_adjoint;
end
img.fwd_model = mdl;
img.elem_data= ones(n_el,1) * conduc ;
Jnp= calc_jacobian(img);

mdl.jacobian = @jacobian_perturb;
Jp1= calc_jacobian(img);

img.elem_data= ones(n_el,1) * conduc ;
Jp2= zeros(size(Jnp));
delta= 1e-8;
for i=1:n_el
   fsol_h= fwd_solve(img);
   img.elem_data(i) = img.elem_data(i) + delta;
   fsol_i= fwd_solve(img);
   Jp2(:,i) = (fsol_i.meas - fsol_h.meas)/delta; 
end
   
mdl.solve = @fwd_solve_1st_order;
mdl.jacobian = @jacobian_adjoint;
Jaa= calc_jacobian(img);

fprintf('Jacobians: Cols by Jaa, Jnp, Jp1, Jp2\n');
JJ = [Jaa;Jnp;Jp1;Jp2];
disp(JJ(:,1:6))
