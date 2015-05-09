function num = num_nodes( mdl );
% NUM_NODES: number of elemnts in a (fwd or inv model or image)
% n_nodes = num_nodes( mdl );

% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

switch mdl.type
  case 'image';      fmdl = mdl.fwd_model;
  case 'inv_model';  fmdl = mdl.fwd_model;
  case 'fwd_model';  fmdl = mdl;
  otherwise;
      error('can''t process model of type %s', mdl.type );
end

num = size(fmdl.nodes,1);

function do_unit_test
   mdl = mk_common_model('a2c2',8);
   ne = num_nodes( mdl );
   unit_test_cmp('test1',ne,41);

   ne = num_nodes( mdl.fwd_model );
   unit_test_cmp('test2',ne,41);

   ne = num_nodes( mk_image( mdl ));
   unit_test_cmp('test3',ne,41);

   mdl = mk_common_model('n3r2',[16,2]);
   ne = num_nodes( mk_image( mdl ));
   unit_test_cmp('test4',ne,252);
