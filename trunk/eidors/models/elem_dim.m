function num = elem_dim( mdl );
% ELEM_DIM: dimension of elements in space (are elements in 2D or 3D space)
%  Note, elem_dim is not the same as the model dimention (can have 2D models in 3D space)
% num_dim = elem_dim( mdl );
%    = 2 if 2D, 3 if 3D

% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

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

num = size(fmdl.elems,2)-1;

function do_unit_test
   mdl = mk_common_model('a2c2',8);
   ne = elem_dim( mdl );
   unit_test_cmp('test1', ne, 2);

   mdl = mk_common_model('n3r2',[16,2]);
   ne = elem_dim( mk_image( mdl ));
   unit_test_cmp('test2', ne, 3);
