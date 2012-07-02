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

if( size(img.elem_data,1) == num_elems(fwd_model))
    %For non-dual models
    elem_data = img.elem_data; 
else
    %For dual models
    if isfield(fwd_model, 'coarse2fine')
        elem_data = fwd_model.coarse2fine * img.elem_data;
    end  

% Fixme: is 'background' the right parameter here?
    %If background image is known
    if isfield(fwd_model, 'background')
        elem_data = elem_data + fwd_model.background;
    end
end

elem_sigma = kron( elem_data(:), mdl_dim(fwd_model) );

ES= ones(lFC,1);
ES(1:length(elem_sigma))= elem_sigma;
ES= spdiags(ES,0,lFC,lFC);

s_mat.E= FC' * ES * FC;

function do_unit_test
   imdl=  mk_common_model('a2c2',16);
   img=  mk_image(imdl);
   S1 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat1', S1.E(1,:), [8,-2,-2,-2,-2,zeros(1,36)],1e-14);

   img.elem_data([1:16]) = 2;
img.calc_colours.clim = 4;
show_fem(img,[0,0,3]);
   S2 = system_mat_1st_order(img.fwd_model,img);
   unit_test_cmp('sys_mat2', S2.E(1,:), 2*[8,-2,-2,-2,-2,zeros(1,36)],1e-14);
   
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
