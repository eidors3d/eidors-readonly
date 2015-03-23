function num = num_elecs( mdl );
% NUM_ELECS: number of electrodes attached to model
% num_elecs = num_elecs( mdl );

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

if isfield(fmdl,'electrode')
   num = length(fmdl.electrode);
else 
   num = 0;
end


function do_unit_test
   mdl = mk_common_model('a2c2',8);
   ne = num_elecs( mdl );
   ok='fail'; if ne==8; ok='ok'; end; fprintf('test1: %10s\n',ok);

   ne = num_elecs( mdl.fwd_model );
   ok='fail'; if ne==8; ok='ok'; end; fprintf('test2: %10s\n',ok);

   ne = num_elecs( mk_image( mdl ));
   ok='fail'; if ne==8; ok='ok'; end; fprintf('test3: %10s\n',ok);

   mdl.fwd_model.electrode = struct([]);
   ne = num_elecs( mdl );
   ok='fail'; if ne==0; ok='ok'; end; fprintf('test4: %10s\n',ok);

   mdl.fwd_model = rmfield(mdl.fwd_model,'electrode');
   ne = num_elecs( mdl );
   ok='fail'; if ne==0; ok='ok'; end; fprintf('test5: %10s\n',ok);

   mdl = mk_common_model('n3r2',[16,2]);
   ne = num_elecs( mk_image( mdl ));
   ok='fail'; if ne==32; ok='ok'; end; fprintf('test6: %10s\n',ok);
