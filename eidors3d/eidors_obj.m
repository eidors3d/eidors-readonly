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
% USAGE: as an accessor
%   gets the version of obj from the eidors_obj database
%     obj  = eidors_obj('get',obj );
%  or
%     obj  = eidors_obj( obj );
%
% USAGE: to log status messages
%     obj  = eidors_obj('msg', message, loglevel )
%
% this will get or set the values of cached properties of the object.
%

% TODO: 
%   1. add code to delete old objects
%   2. accessors and setters of the form 'prop1.subprop1'

% FIXME: this variable must be global, rather than persistent
% because the broken Matlab syntax does not allow persistent
% variables to be passed to subfunctions
% global eidors_objects
% (Short circuit boolean removed for compatibility with Matlab 6.1 (R12.1) WRBL 22/02/2004)
% Converted eidors_objects.(x) to getfield or setfield WRBL 22/02/2004
if nargin==0 | ~isstr(type)
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

%   obj = eidors_objects.( obj_id );
    obj = getfield(eidors_objects, obj_id )

function obj = set_obj( obj, varargin );
   global eidors_objects
   obj_id= test_exist( obj );

%  eidors_objects.( obj_id ) = obj;
%  eidors_objects.( obj_id ).cache= []; %clear cache
%  
%  for idx= 1:2:nargin-1
%     eidors_objects.( obj_id ).( varargin{idx} )= ...
%                                 varargin{idx+1};
%  end
%  obj = eidors_objects.( obj_id );

% for octave and matlab 6.1 compatibility
   obj.cache= [];
   for idx= 1:2:nargin-1
      eval(sprintf('obj.%s=varargin{%d};', varargin{idx},idx+1 ));
   end
   eval(sprintf('eidors_objects.%s=obj;',obj_id));

function val= cache_obj( obj, prop, value );
   global eidors_objects
   obj_id= test_exist( obj );

   if nargin==3     %set cache
%     eidors_objects.( obj_id ).cache.( prop ) = value;
      eval(sprintf('eidors_objects.%s.cache.%s=value;', obj_id, prop));
      val= ''; % to satisfy matlab errors
   elseif nargin==2 %get from cache
      try
%       val= eidors_objects.( obj_id ).cache.( prop );
        val = eval(sprintf('eidors_objects.%s.cache.%s;',obj_id,prop));
      catch
         val= [];
      end
   else
      error('nargin wrong calling cache_obj');
   end

function obj= new_obj( type, name, varargin );
   global eidors_objects
   try
      obj_id =  sprintf('obj_%07d', eidors_objects.idnum);
   catch
      eidors_objects.idnum= 10000001;
      obj_id= 'obj_10000001';
   end

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
   obj= set_obj(obj, varargin{:} );

% This function hashes the value of the variable var
% to create an obj_id. The goal is to allow proper caching
% of calculated matrices, by detecting when a previous
% calculation with same parameters has been made
%
% Matlab does not offer any way to create hashes of variables.
% It does not even offer a good way to get the size of the total
% variable storage.  Thus, I choose to save the variable, calculate
% the sha1 hash with sha1sum, and delete the file. Given a modern
% OS, file data should only be cached in memory. The largest overhead
% will be for the instantiation of the process
function obj_id= obj_hash( var )
   global eidors_objects;

   tmpnam= [tempname ,'.mat'];
   save('-v6',tmpnam,var);
   delete(tmpnam);


