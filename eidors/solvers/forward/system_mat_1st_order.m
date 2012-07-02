function s_mat= system_mat_1st_order( fwd_model, img)
% SYSTEM_MAT_1ST_ORDER: SS= system_mat_1st_order( fwd_model, img)
% Calc system matrix for Andy Adler's EIT code
% fwd_model = forward model
% img       = image background for system matrix calc
% s_mat.E = CC' * SS * conductivites * CC;
% where:
%   SS  = Unconnected system Matrix
%   CC  = Connectivity Matrix

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return; end

FC= system_mat_fields( fwd_model);
lFC= size(FC,1);

elem_data = check_elem_data(fwd_model, img);
elem_sigma = kron( elem_data, ones(elem_dim(fwd_model),1) );
elem_sigma(end+1:lFC) = 1; % add ones for CEM

ES= spdiags(elem_sigma,0,lFC,lFC);

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
       otherwise; error(['system_mat_1st_order: provided elem_data ' ...
            ' (sz=%d) does not match c2f (sz=%d %d)'], sz_elem_data(1), sz_c2f);
     end
   else
     if sz_elem_data(1) ~= num_elems(fwd_model)
       error(['system_mat_1st_order: provided elem_data (sz=%d) does ' ...
          ' not match fwd_model (sz=%d)'], sz_elem_data(1), num_elems(sz_c2f));
     end
   end


% FIXME: is 'background' the right parameter here?
% There is some need to provide a background, especially
%  for absolute solutions with a smaller c2f model, 
%  but this doesn't seem the best way. It's here for now
%  (as of July, 2012), but we plan to develop a better
%  solution.
    if isfield(fwd_model, 'background')
        elem_data = elem_data + fwd_model.background;
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
      S2 = system_mat_1st_order(img.fwd_model,img);
      unit_test_cmp('sys_mat: error', 1,0);
   catch
      unit_test_cmp('sys_mat: error', 1,1);
   end

   cmdl = mk_circ_tank(4,[],0);
   c2f = mk_coarse_fine_mapping( img.fwd_model, cmdl); % spy(c2f)
   img.fwd_model.coarse2fine = c2f;
   img.elem_data = ones(size(cmdl.elems,1),1);
   S3 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('c2f #1', S1.E, S3.E,1e-14);

   cmdl = mk_circ_tank(2,[],0);
   c2f = mk_coarse_fine_mapping( img.fwd_model, cmdl); % spy(c2f)
   img.fwd_model.coarse2fine = c2f;
   img.elem_data = ones(size(cmdl.elems,1),1);
   S4 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('c2f #2', S1.E, S4.E,1e-14);

