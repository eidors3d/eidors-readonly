function num = mdl_dim( mdl );
% MDL_DIM: dimension of model space (are nodes in 2D or 3D space)
%  Note, mdl_dims is not the same as the element dimentions
% num_dim = mdl_dim( mdl );
%    = 2 if 2D, 3 if 3D

% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

if ~isfield(mdl,'type') && isfield(mdl,'nodes') % Set default for this case
   mdl.type = 'fwd_model';
end

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
   ne = mdl_dim( mdl );
   ok='fail'; if ne==2; ok='ok'; end; fprintf('test1: %10s\n',ok);

   mdl = mk_common_model('n3r2',16);
   ne = mdl_dim( mk_image( mdl ));
   ok='fail'; if ne==3; ok='ok'; end; fprintf('test2: %10s\n',ok);
