function s_mat= system_mat_1st_order( fwd_model, img)
% SYSTEM_MAT_1ST_ORDER: SS= system_mat_1st_order( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% img.fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

if nargin == 1
   img= fwd_model;
elseif  strcmp(getfield(warning('query','EIDORS:DeprecatedInterface'),'state'),'on')
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling SYSTEM_MAT_1ST_ORDER with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
fwd_model= img.fwd_model;

% check physics
if isfield(img,'current_physics') && ~isempty(img.current_physics) ... 
        && ~strcmp(img.current_physics,'conductivity')
    error('system_mat_1st_order does not work for %s',img.current_physics);
end

FC= system_mat_fields( fwd_model);
lFC= size(FC,1);
ED = elem_dim(fwd_model);
lNE= ED*num_elems(fwd_model);

elem_data = check_elem_data(fwd_model, img);
if size(elem_data,3) == 1
% Scalar conductivity == isotropic
   elem_sigma = kron( elem_data, ones(ED,1) );
   elem_sigma(end+1:lFC) = 1; % add ones for CEM

   ES= spdiags(elem_sigma,0,lFC,lFC);
else
   switch elem_dim(fwd_model)
     case 2;
       idx = 1:2:lNE;
       ES= sparse([idx,idx+1,idx,idx+1]', ...
                  [idx,idx,idx+1,idx+1]', elem_data(:), lFC,lFC);
       
       ES(lNE+1:lFC,lNE+1:lFC) = speye(lFC-lNE);
     case 3;
       idx = 1:3:lNE;
       ES= sparse([idx,idx+1,idx+2,idx,idx+1,idx+2,idx,idx+1,idx+2]', ...
                  [idx,idx,idx,idx+1,idx+1,idx+1,idx+2,idx+2,idx+2]', ...
                   elem_data(:), lFC,lFC);
       
       ES(lNE+1:lFC,lNE+1:lFC) = speye(lFC-lNE);
     otherwise; 
       error('%d D anisotropic elements not implemented', elem_dim(fwd_model));
   end

end

s_mat.E= FC' * ES * FC;

function elem_data = check_elem_data(fwd_model, img);
   elem_data = img.elem_data; 
   sz_elem_data = size(elem_data);
   if sz_elem_data(2) ~= 1;
      error('system_mat_1st_order: can only solve one image (sz_elem_data=%)', ...
            sz_elem_data);
   end

   if isfield(fwd_model, 'coarse2fine');
     c2f = fwd_model.coarse2fine;
     sz_c2f = size(c2f);
     switch sz_elem_data(1)
       case sz_c2f(1); % Ok     
       case sz_c2f(2); elem_data = c2f * elem_data;
          if isfield(fwd_model, 'background')
              elem_data = elem_data + fwd_model.background;
          end

       otherwise; error(['system_mat_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model)
       error(['system_mat_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end


% Need tests for each error case
% Need tests for background

function do_unit_test
   imdl=  mk_common_model('a2c2',16);
   img=  mk_image(imdl);
   S1 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat1', S1.E(1,:), [4,-1,-1,-1,-1,zeros(1,36)],1e-14);

   img.elem_data([1:16]) = 2;
   img.calc_colours.clim = 4; show_fem(img,[0,0,3]);

   S2 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat2', S2.E(1,:), 2*[4,-1,-1,-1,-1,zeros(1,36)],1e-14);

   imdl=  mk_common_model('a2C2',16);
   img=  mk_image(imdl);
   S1 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat1', S1.E(1,:), [4,-1,-1,-1,-1,zeros(1,52)],1e-14);

   img.elem_data([1:16]) = 2;
   img.calc_colours.clim = 4; show_fem(img,[0,0,3]);

   S2 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat2', S2.E(1,:), 2*[4,-1,-1,-1,-1,zeros(1,52)],1e-14);
   
   idx = 41-(0:15);
   unit_test_cmp('sys_mat4', S2.E(idx,idx),   S1.E(idx,idx),1e-14);
   idx = 1:5;
   unit_test_cmp('sys_mat3', S2.E(idx,idx), 2*S1.E(idx,idx),1e-14);

   img.elem_data([1:36]) = 2;
   S2 = system_mat_1st_order(img.fwd_model,img);
   idx = 1:5;
   unit_test_cmp('sys_mat3', S2.E(idx,idx), 2*S1.E(idx,idx),1e-14);

   img.elem_data = ones(63,1);
   try
     % this should give an error, since elem_data is wrong size
      S2 = system_mat_1st_order(img.fwd_model,img);
      unit_test_cmp('sys_mat: test for size error', 1,0);
   catch
      unit_test_cmp('sys_mat: test for size error', 1,1);
   end

   imdl=  mk_common_model('a2c2',16);
   cmdl = mk_circ_tank(4,[],0);
   c2f = mk_coarse_fine_mapping( img.fwd_model, cmdl); % spy(c2f)
   img.fwd_model.coarse2fine = c2f;
   img.elem_data = ones(size(cmdl.elems,1),1);
   S3 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('c2f #1', S1.E, S3.E,1e-14);

   cmdl = mk_circ_tank(2,[],0);
   c2f = mk_coarse_fine_mapping( img.fwd_model, cmdl); % spy(c2f)
   c2f = c2f./(sum(c2f,2) * ones(1,size(c2f,2))); % FIX
   img.fwd_model.coarse2fine = c2f;
   img.elem_data = ones(size(cmdl.elems,1),1);
   S4 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('c2f #2', S1.E, S4.E,1e-14);

   test_2d_resistor
   test_3d_resistor

function test_2d_resistor
   nn= 12;     % number of nodes
   ww=2;       % width = 4
   conduc=  .4;% conductivity in Ohm-meters
   current= 4;  % Amps
   z_contact= 9e-1;
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
   mdl.normalize_measurements = 0;
   %show_fem(mdl);

   % analytical solution
   wid_len= max(mdl.nodes) - min(mdl.nodes);
   R = wid_len(2) / wid_len(1) / conduc + 2*z_contact/scale;

   V= current*R;
   %fprintf('Solver %s: %f\n', 'analytic', V);

   % AA_SOLVER
   mdl.solve = @fwd_solve_1st_order;
   mdl.system_mat = @system_mat_1st_order;
   img = mk_image(mdl, conduc);
   fsol= fwd_solve(img);
   %fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('test 2d vs analytic', fsol.meas, V, 1e-11);

   % NP_SOLVER
   mdl.solve = @np_fwd_solve;
   mdl.system_mat = @np_calc_system_mat;
   img = mk_image(mdl, conduc);
   fsol= fwd_solve(img);
   %fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('np_fwd_solve vs analytic', fsol.meas, V, 1e-11);

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

   mdl.solve = @fwd_solve_1st_order;
   mdl.system_mat = @system_mat_1st_order;

   img= mk_image(mdl, conduc);

   fsol= fwd_solve(img);
%  fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('test 2d vs analytic', fsol.meas, V, 1e-11);

   mdl.solve = @np_fwd_solve;
   mdl.system_mat = @np_calc_system_mat;
   mdl.misc.perm_sym= '{n}';
   img= mk_image(mdl, conduc);

   fsol= fwd_solve(img);
%  fprintf('Solver %s: %f\n', fsol.name, fsol.meas);
   unit_test_cmp('np_fwd_solve vs analytic', fsol.meas, V, 1e-11);
