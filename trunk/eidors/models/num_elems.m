function num = num_elems( mdl );
% NUM_ELEMS: number of elemnts in a (fwd or inv model or image)
% num_elems = num_elems( mdl );

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

num = size(fmdl.elems,1);


function do_unit_test
   mdl = mk_common_model('a2c2',8);
   ne = num_elems( mdl );
   unit_test_cmp('test1',ne, 64);

   ne = num_elems( mdl.fwd_model );
   unit_test_cmp('test2',ne, 64);

   ne = num_elems( mk_image( mdl ));
   unit_test_cmp('test3',ne, 64);

   mdl = mk_common_model('n3r2',[16,2]);
   ne = num_elems( mk_image( mdl ));
   unit_test_cmp('test4',ne, 828);
