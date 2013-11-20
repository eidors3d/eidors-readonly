function num = ls_elem_dim( mdl );

if ~isfield(mdl,'type') && isfield(mdl,'elems') % Set default for this case
   mdl.type = 'fwd_model';
end

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

n = size(fmdl.elems,2);
if n==4
    num=2;
elseif n==8
    num=3;
end


