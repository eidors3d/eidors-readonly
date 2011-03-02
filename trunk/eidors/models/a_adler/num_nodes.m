function num = n_nodes( mdl );
% N_NODES: number of elemnts in a (fwd or inv model or image)
% num_nodes = n_nodes( mdl );

% $Id$

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

num = size(fmdl.nodes,1);
