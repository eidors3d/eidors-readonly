function img = physics_param_mapper(img, rev)

if nargin < 2
   rev  = 0;
end
if rev
   img = reverse_mapping(img);
else
   img = forward_mapping(img);
end
   
function img = forward_mapping(img)
   if ~isfield(img, 'physics_param_mapper')
      jnk = physics_data_mapper(img);
      try 
         img.params = jnk.node_data;
         img.current_params = 'node_data';
      catch
         img.params = jnk.elem_data;
         img.current_params = 'elem_data';
      end
   elseif iscell(img.physics_param_mapper)
      params = [];
      for i = 1:length(img.physics_param_mapper)
         params = [params; eval(['img.' img.physics_param_mapper{i}])];
      end
      img.params = params;
   end
      
   
function img = reverse_mapping(img)
   if ~isfield(img, 'physics_param_mapper')
      if ~isfield(img, 'current_params')
         error('Need img.current_params');
      end
      switch img.current_params
         case 'node_data'
            img.node_data = img.params;
         case 'elem_data'
            img.elem_data = img.params;
         otherwise
            error('huh?');
      end
      img = physics_data_mapper(rmfield(img,'params'));
      img.current_params = [];
   elseif iscell(img.physics_param_mapper)
      start = 1;
      for i = 1:length(img.physics_param_mapper)
         len = size(eval(['img.' img.physics_param_mapper{i}]),1);
         stop = start+len-1;
         eval(['img.' img.physics_param_mapper{i} ' = img.params(start:stop,:);']);
         start = len + 1;
      end
  end
      