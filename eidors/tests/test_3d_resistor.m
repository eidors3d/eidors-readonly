% Create 3D model of a Rectangular resistor
% $Id$

ll=5; % length
ww=1; % width
hh=1; % height
conduc= .13;  % conductivity in Ohm-meters
current= 4;  % Amps
z_contact= 1e-2;
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

mdl.solve = @fwd_solve_1st_order;
mdl.system_mat = @system_mat_1st_order;

n_el = size(mdl.elems,1);
img= eidors_obj('image','3D rectangle', ...
      'elem_data', ones(n_el,1) * conduc, ...
      'fwd_model', mdl); 

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

mdl.solve = @np_fwd_solve;
mdl.system_mat = @np_calc_system_mat;
mdl.misc.perm_sym= '{n}';
mdl.normalize_measurements = 0;
img= eidors_obj('image','3D rectangle', ...
      'elem_data', ones(size(mdl.elems,1),1) * conduc, ...
      'fwd_model', mdl); 

fsol= fwd_solve(img);
fprintf('Solver %s: %f\n', fsol.name, fsol.meas);

% analytical solution
R = ll / ww / hh / scale/ conduc + 2*z_contact/scale^2;

V= current*R;
fprintf('Solver %s: %f\n', 'analytic', V);


% NOW CALCULATE THE ANALYTICAL JACOBIAN
mdl.solve = @np_fwd_solve;
mdl.jacobian = @np_calc_jacobian;
img.fwd_model = mdl;
img.elem_data= ones(n_el,1) * conduc ;
Jnp= calc_jacobian(img);

mdl.jacobian = @jacobian_perturb;
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
   
mdl.solve = @fwd_solve_1st_order;
mdl.jacobian = @jacobian_adjoint;
Jaa= calc_jacobian(mdl,img);

[Jaa;Jnp;Jp1;Jp2]
