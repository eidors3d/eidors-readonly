function obj_id= eidors_obj(type,name, varargin );
% EIDORS_OBJ: 'constructor' to create a eidors structure
% USAGE: as a constructor
%     obj  = eidors_obj(type,name,prop1,value1, prop1, value2, ...)
%
%     type:  obj type: fwd_model, inv_model, data, image
%     name:  text string identifier for model (may be '')
%     prop1, value1, etc: properites and values for object
%
%    example: fwd_mdl = ...
%             eidors_obj('fwd_model','My FWD MODEL', ...
%                        'nodes', NODES, 'elems', ELEMS, ...
%
% OR construct from structure
%     obj  = eidors_obj(type,obj);
%
%     example: fwd_mdl.nodes = NODES; .... %etc
%              fwd_mdl = eidors_obj('fwd_model',fwd_mdl);
%
% USAGE: to set values
%     obj  = eidors_obj('set',obj,prop1,value1, prop1, value2, ...)
%
% this will set the values of properties of the object. At the
% same time, any cached values will be erased (because they may
% depend on older properties of the model)
%   
%    example:
%             eidors_obj('set',fwd_mdl, 'nodes', NEW_NODES);
%
% USAGE: to cache values
%     obj  = eidors_obj('cache',obj, cacheprop1,value1, ...)
%
% this will get or set the values of cached properties of the object.
%
%    example: % set jacobian
%             eidors_obj('cache',fwd_mdl, 'jacobian', J):
%
%    example: % get jacobian or '[]' if not set
%             J= eidors_obj('cache',fwd_mdl, 'jacobian'):
%
%
% USAGE: as an accessor
%     obj  = eidors_obj('get',objid, prop1);
%  or
%     obj  = eidors_obj('get',objid );
%  or
%     obj  = eidors_obj( objid );

% TODO: 
%   1. add code to delete old objects
%   2. accessors and setters of the form 'prop1.subprop1'

% FIXME: this variable must be global, rather than persistent
% because the broken Matlab syntax does not allow persistent
% variables to be passed to subfunctions
% global eidors_objects

if nargin==0 || ~isstr(type)
   error('cannot call eidors_obj with no arguments');
end

if nargin==1
   obj_id= get_obj( name );
elseif nargin>1
   if     strcmp( type, 'get')
      obj_id= get_obj( name, varargin{:} );
   elseif strcmp( type, 'set')
      obj_id= set_obj( name, varargin{:} );
   elseif strcmp( type, 'cache')
      obj_id= cache_obj( name, varargin{:} );
   else
      obj_id= new_obj( type, name, varargin{:} );
   end
end

function obj_id= test_exist( obj );
   global eidors_objects
   obj_id = obj.id;   
   if ~isfield( eidors_objects , obj_id)
      error(sprintf('EIDORS object %s does not exist', obj_id));
   end

function obj = get_obj( obj, varargin );
   global eidors_objects
   obj_id= test_exist( obj );

   if nargin==1
      obj = eidors_objects.( obj_id );
   elseif nargin==2
      obj = eidors_objects.( obj_id ).( varargin{1} );
   else
      error('get_obj cannot interpret input');
   end

function obj = set_obj( obj, varargin );
   global eidors_objects
   obj_id= test_exist( obj );

   eidors_objects.( obj_id ) = obj;
   eidors_objects.( obj_id ).cache= []; %clear cache
   for idx= 1:2:nargin-1
      eidors_objects.( obj_id ).( varargin{idx} )= ...
                                  varargin{idx+1};
   end

   obj = eidors_objects.( obj_id );
     
function val= cache_obj( obj, prop, value );
   global eidors_objects
   obj_id= test_exist( obj );

   if nargin==3     %set cache
      eidors_objects.( obj_id ).cache.( prop ) = value;
      val= ''; % to satisfy matlab errors
   elseif nargin==2 %get from cache
      try
         val= eidors_objects.( obj_id ).cache.( prop );
      catch
         val= [];
      end
   else
      error('nargin wrong calling cache_obj');
   end

function obj= new_obj( type, name, varargin );
   global eidors_objects
   if isempty(eidors_objects)
       eidors_objects.idnum = 10000001;
   end
   obj_id =  sprintf('obj_%07d', eidors_objects.idnum);
   eidors_objects.idnum  = eidors_objects.idnum + 1;

   % if called with obj, work as a copy-constructor
   if isstruct(name)
      obj= name;
      try
         name= obj.name;
      catch
         name= 'unknown';
      end
   end

   obj.type = type;
   obj.id   = obj_id;
   obj.name = name;

   eval(sprintf('eidors_objects.%s=obj;', obj_id)); % create in store
   obj.name = name;
   obj.id   = obj_id;
   obj= set_obj(obj, varargin{:} );
