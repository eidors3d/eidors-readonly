function num = mdl_dims( mdl );
% MDL_DIMS: dimension of model space (are nodes in 2D or 3D space)
%  Note, mdl_dims is not the same as the element dimentions
% num_dims = n_dims( mdl );
% num_dims = 2 if 2D, 3 if 3D

% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

num = size(fmdl.nodes,2);

function do_unit_test
   mdl = mk_common_model('a2c2',8);
   ne = mdl_dims( mdl );
   ok='fail'; if ne==2; ok='ok'; end; fprintf('test1: %10s\n',ok);

   mdl = mk_common_model('n3r2',16);
   ne = mdl_dims( mk_image( mdl ));
   ok='fail'; if ne==3; ok='ok'; end; fprintf('test2: %10s\n',ok);
