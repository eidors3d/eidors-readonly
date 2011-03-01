function num = n_dims( mdl );
% N_DIMS: dimension of FEM in a (fwd or inv model or image)
% num_dims = n_dims( mdl );
% num_dims = 2 if 2D, 3 if 3D

% $Id$

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

num = size(fmdl.nodes,2);
