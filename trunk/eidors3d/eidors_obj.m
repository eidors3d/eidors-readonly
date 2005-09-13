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
%          eidors_obj('set-cache',obj, cacheprop1,value1, dep_objs, ...)
%     obj= eidors_obj('get-cache',obj, cacheprop1, dep_objs, ...)
%
% this will get or set the values of cached properties of the object.
%
%    example: % set jacobian
%             eidors_obj('set-cache',fwd_mdl, 'jacobian', J):
%
%    example: % get jacobian or '[]' if not set
%             J= eidors_obj('get-cache',fwd_mdl, 'jacobian'):
%
% However, in some cases, such as the Jacobian, the value depends
% on other objects, such as the image background. In this case, use
%
%    example: % set jacobian
%             eidors_obj('set-cache',fwd_mdl, 'jacobian', J, homg_img):
%
%    example: % get jacobian or '[]' if not set
%             J= eidors_obj('get-cache',fwd_mdl, 'jacobian', homg_img):
%
%
% USAGE: to log status messages
%     obj  = eidors_obj('msg', message, loglevel )
%
% this will get or set the values of cached properties of the object.
%

% $Id: eidors_obj.m,v 1.19 2005-09-13 02:27:29 aadler Exp $
% TODO: 
%   1. add code to delete old objects
%   2. accessors and setters of the form 'prop1.subprop1'

% FIXME: this variable must be global, rather than persistent
% because the broken Matlab syntax does not allow persistent
% variables to be passed to subfunctions
% global eidors_objects
% (Short circuit boolean removed for compatibility with Matlab 6.1 (R12.1) WRBL 22/02/2004)
% Converted eidors_objects.(x) to getfield or setfield WRBL 22/02/2004

% FIXME!!! Bloody Matlab doesn't save its variables consistently
% Testcode:
%   for i=1:20; t1=mk_stim_patterns(n_elec, n_rings, '{ad}','{ad}',options, 10); save(sprintf('t%02d.mat',i),'t1');end
% for i in t*.mat ; do dd if=$i  skip=100 ibs=1 | xxd > $i.xxd ; done
% for i in *.xxd ; do sha1sum $i ; done
% 6! different saved formats exist. Damn

if nargin==0 | ~isstr(type)
   error('cannot call eidors_obj with no arguments');
end

if nargin==1
   obj_id= get_obj( name );
elseif nargin>1
   if     strcmp( type, 'set')
      obj_id= set_obj( name, varargin{:} );
   elseif strcmp( type, 'get-cache')
      obj_id= get_cache_obj( name, varargin{:} );
   elseif strcmp( type, 'set-cache')
      set_cache_obj( name, varargin{:} );
      obj_id= []; % quiet matlab errors
   else
      obj_id= new_obj( type, name, varargin{:} );
   end
end

function obj = set_obj( obj, varargin );
   global eidors_objects
% we choose not to cache data and images because this will
% tend to fill up the workspace.
% TODO: make the DONT_CACHE list use configuable
   DONT_CACHE= {'data','image'};

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

   if any(strcmp( obj.type, DONT_CACHE)); return; end 
   obj_id= calc_obj_id( obj );
      
% set the obj_id into the obj after calculating the hash value
   obj.id= obj_id;
   if ~isfield(eidors_objects,obj_id)
      eval(sprintf('eidors_objects.%s=obj;',obj_id));
   end

function val= get_cache_obj( obj, prop, dep_obj1 );
   global eidors_objects
   val= [];

   try
      obj_id= obj.id;
   catch
      return; % Can't cache onto non-existant objects
   end

   if nargin==3
      try
          prop= [prop,'_', dep_obj1.id];
      catch
          prop= [prop,'_', calc_obj_id(dep_obj1)];
      end
   end

   try
%     val= eidors_objects.( obj_id ).cache.( prop );
      val = eval(sprintf('eidors_objects.%s.cache.%s;',obj_id,prop));
   end

function set_cache_obj( obj, prop, value, dep_obj1 );
   global eidors_objects

   try
      obj_id= obj.id;
   catch
      return; % Can't cache onto non-existant objects
   end

   if nargin==4
      try
          prop= [prop,'_', dep_obj1.id];
      catch
          prop= [prop,'_', calc_obj_id(dep_obj1)];
      end
   end


%  eidors_objects.( obj_id ).cache.( prop ) = value;
   eval(sprintf('eidors_objects.%s.cache.%s=value;', obj_id, prop));

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
   if ~isfield(eidors_objects,'hash_type')
       test_for_hashtypes;
   end

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
 
   if eidors_objects.hash_type==1; 
       obj_id= matrix_id(tmpnam);
   elseif eidors_objects.hash_type > 1e6
%if hashing code is unavailable, then disable caching function
       obj_id= sprintf('id_%040d', eidors_objects.hash_type );
       eidors_objects.hash_type= eidors_objects.hash_type + 1;
   else
       error('hash_type value unrecognized');
   end
   delete(tmpnam);
%  fprintf('creating [%s] object %s\n',var.type,obj_id);

% Since we use a dynamically loaded function to test
% for hashtypes, all sorts of things can go wrong
%   - The compilers or libraries may not be available
%
% We test for:
%   1. Existance of 'matrix_id' mex file
%   2. Existance of utilities 'sha1sum' to do this
% 
% If nothing exists, then the best we can do is to
% disable hashing completely. We set hashtype to a
% number and increment it each time   
function test_for_hashtypes
   global eidors_objects; 
   if exist('matrix_id')==3 % MEX-file on MATLAB's search path
      eidors_objects.hash_type = 1;
   else
      eidors_objects.hash_type= 1e6+1; 
   end

