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
% USAGE: to log status messages
%     obj  = eidors_obj('msg', message, loglevel )
%
% this will get or set the values of cached properties of the object.
%

% $Id: eidors_obj.m,v 1.13 2005-02-23 14:43:48 aadler Exp $
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
   if     strcmp( type, 'set')
      obj_id= set_obj( name, varargin{:} );
   elseif strcmp( type, 'cache')
      obj_id= cache_obj( name, varargin{:} );
   else
      obj_id= new_obj( type, name, varargin{:} );
   end
end

function obj = set_obj( obj, varargin );
   global eidors_objects

   if strcmp(obj.type, 'inv_model'); keyboard; end
   try
      old_obj_id= obj.id;
      % If we're modifying an old object, then remove the old version
      eidors_objects= rmfield(eidors_objects, old_obj_id);
   end
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
   obj_id= calc_obj_id( obj );
      
% set the obj_id into the obj after calculating the hash value
   obj.id= obj_id;
   if ~isfield(eidors_objects,obj_id)
      eval(sprintf('eidors_objects.%s=obj;',obj_id));
   end

function val= cache_obj( obj, prop, value );
   global eidors_objects

   try
      obj_id= obj.id;
   catch
      error('EIDORS object does not exist');
   end


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

   if isstruct(name)
      obj= name;
      try
         name= obj.name;
      catch
         name= 'unknown';
      end
   end

   obj.type = type;
   obj.name = name;
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
function obj_id= calc_obj_id( var )
   global eidors_objects;

% the obj_id does not depend on the cache or id
   try; var = rmfield(var,'cache'); end
   try; var = rmfield(var,'id');    end

   tmpnam= [tempname ,'.mat'];

% This is another awful example of matlab's
% versionitis. Why can't they stay compatible?
   if version(1)>='7' 
       save(tmpnam,'var','-v6');
   elseif exist('OCTAVE_VERSION')
       save(tmpnam,'-mat','var');
   else
       save(tmpnam,'var');
   end

   obj_id= matrix_id(tmpnam);
   delete(tmpnam);
%  fprintf('creating [%s] object %s\n',var.type,obj_id);


