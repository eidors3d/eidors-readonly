function num = n_elems( mdl );
% N_ELEMS: number of elemnts in a (fwd or inv model or image)
% num_elems = n_elems( mdl );

% $Id$

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

num = size(fmdl.elems,1);
