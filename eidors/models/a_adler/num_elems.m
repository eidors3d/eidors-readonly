function num = num_elems( mdl );
% NUM_ELEMS: number of elemnts in a (fwd or inv model or image)
% num_elems = num_elems( mdl );

% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

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
   ok='fail'; if ne==64; ok='ok'; end; fprintf('test1: %10s\n',ok);

   ne = num_elems( mdl.fwd_model );
   ok='fail'; if ne==64; ok='ok'; end; fprintf('test2: %10s\n',ok);

   ne = num_elems( mk_image( mdl ));
   ok='fail'; if ne==64; ok='ok'; end; fprintf('test3: %10s\n',ok);

   mdl = mk_common_model('n3r2',16);
   ne = num_elems( mk_image( mdl ));
   ok='fail'; if ne==828; ok='ok'; end; fprintf('test4: %10s\n',ok);
