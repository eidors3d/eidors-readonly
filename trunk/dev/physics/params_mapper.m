function img = params_mapper(img, rev)

if nargin < 2
   rev  = 0;
end
if rev
   img = reverse_mapping(img);
else
   img = forward_mapping(img);
end
   
function img = forward_mapping(img)
   if ~isfield(img, 'params_mapper')
      jnk = data_mapper(img);
      phys = [];
%       if ~isempty(jnk.current_params)
%         phys = [jnk.current_params '.'];
%       end
      try 
         img.inv_params = jnk.node_data;
         img.current_inv_params = [phys 'node_data'];
%          img = rmfield(img, 'node_data');
      catch
         img.inv_params = jnk.elem_data;
         img.current_inv_params = [phys 'elem_data'];
%          img = rmfield(img, 'elem_data');
      end
      img.current_params = jnk.current_params; % needed to reverse
   elseif iscell(img.params_mapper)
      inv_params = [];
      for i = 1:length(img.params_mapper)
         inv_params = [inv_params; eval(['img.' img.params_mapper{i}])];
      end
      img.inv_params = inv_params;
      img.current_inv_params = img.params_mapper;
   elseif isa(img.params_mapper,'function_handle')
      img.inv_params = feval(img.params_mapper,img);
      img.current_inv_params = fun2str(img.params_mapper);
   end
      
function img = reverse_mapping(img)
   if ~isfield(img, 'params_mapper')
      if ~isfield(img, 'current_inv_params')
         error('Need img.current_inv_params');
      end
      switch img.current_inv_params
         case 'node_data'
            img.node_data = img.inv_params;
         case 'elem_data'
            img.elem_data = img.inv_params;
         otherwise
            error('huh?');
      end
      img = data_mapper(rmfield(img,'inv_params'),1);
      img.current_inv_params = [];
   elseif iscell(img.params_mapper)
      start = 1;
      for i = 1:length(img.params_mapper)
         len = size(eval(['img.' img.params_mapper{i}]),1);
         stop = start+len-1;
         try
            eval(['img.' img.params_mapper{i} ' = img.inv_params(start:stop,:);']);
         catch
            stop = start + size(img.fwd_model.coarse2fine,2)-1;
            eval(['img.' img.params_mapper{i} ' = img.inv_params(start:stop,:);']);
         end
         start = stop + 1;
      end
      img = rmfield(img,'inv_params'); % remove parameters after mapping
      img.current_inv_params = [];
   elseif isa(img.params_mapper,'function_handle')
      img.inv_params = feval(img.params_mapper,img,1);
  end
