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
%        eidors_obj('fwd_model','My FWD MODEL', ...
%                    'nodes', NODES, 'elems', ELEMS, ...
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
%        eidors_obj('set',fwd_mdl, 'nodes', NEW_NODES);
%
% USAGE: to cache as required
%     obj= eidors_obj('calc-or-cache', obj, funcname, dep_objs ...
%      
% return obj from cache, or calculate it and store in the cache  
%
% USAGE: to cache values
%          eidors_obj('set-cache',obj, cachename,value1, dep_objs, ...)
%     obj= eidors_obj('get-cache',obj, cachename, dep_objs, ...)
%
% this will get or set the values of cached properties of the object.
%
%    example: % set jacobian
%        eidors_obj('set-cache',fwd_mdl, 'jacobian', J):
%
%    example: % get jacobian or '[]' if not set
%        J= eidors_obj('get-cache',fwd_mdl, 'jacobian'):
%
% However, in some cases, such as the Jacobian, the value depends
% on other objects, such as the image background. In this case, use
%
%    example: % set jacobian
%        eidors_obj('set-cache',fwd_mdl, 'jacobian', J, homg_img):
%
%    example: % get jacobian or '[]' if not set
%        J= eidors_obj('get-cache',fwd_mdl, 'jacobian', homg_img):
 

% OLD FUNCTION - Is this still needed
% In many cases, data can be labelled, and explicitly cached to a file:
% To specify the file, use the last parameter of a 'get-cache',
%   'set-cache' or 'calc-or-cache' command
%
%     example: % set jacobian for homg_img
%        eidors_obj('set-cache',fwd_mdl,'jacobian',J, h_img,'homg'): 
%        J= eidors_obj('get-cache',fwd_mdl,'jacobian',h_img,'homg'):
%
%     this will save to 'cachedir/homg.mat' where cachedir is
%        eidors_object.cachedir. Typically, a subdirectory is 
%        created for each algorithm's data
%
% 

% $Id: eidors_obj.m,v 1.33 2005-10-18 18:56:23 aadler Exp $

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
   elseif strcmp( type, 'get-cache')
      obj_id= get_cache_obj( name, varargin{:} );
   elseif strcmp( type, 'set-cache')
      set_cache_obj( name, varargin{:} );
      obj_id= []; % quiet matlab errors
   elseif strcmp( type, 'calc-or-cache')
      val = calc_or_cache(name, varargin{:} );
      obj_id= val;
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
      % unless it contains cached data
      
      obj= rmfield(obj,'id');
      if ~isfield( eidors_objects.old_obj_id, 'cache' )
         eidors_objects= rmfield(eidors_objects, old_obj_id);
      end
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
   for idx= 1:2:nargin-1
      eval(sprintf('obj.%s=varargin{%d};', varargin{idx},idx+1 ));

   end
   if any(strcmp( obj.type, DONT_CACHE)); return; end 

   % to ensure that function handles become part of eidors_var_id
   % FIX this when we update mex file
   
   try 
      if strcmp( class( obj.image_prior.func ), 'function_handle')
          eval(sprintf('obj.misc.function_handles.ipf=''%s'';', ...
                func2str( obj.image_prior.func) ) );
      end
   end

   obj_id= calc_obj_id( obj );
      
   obj.id= obj_id;

% set the obj_id into the obj after calculating the hash value
   if ~isfield(eidors_objects,obj_id)
      eval(sprintf('eidors_objects.%s=obj;',obj_id));
   end

function val= calc_or_cache( obj, funcname, varargin )

   if strcmp( class(funcname), 'function_handle')
      funcstr= func2str(funcname);
   else
      funcstr= funcname;
   end
   
   val = get_cache_obj( obj, funcname, varargin{:} );
   if ~isempty( val );
      eidors_msg([funcstr, ': using cached value'], 2);
   else
      [objlist, cachename]= proc_obj_list( varargin{:} );
      val = feval( funcname, obj, objlist{:} );
      set_cache_obj( obj, funcstr, val, varargin{:} );
      eidors_msg([funcstr, ': setting cached value'], 2);
   end

% val= get_cache_obj( obj, prop, dep_obj1, dep_obj2, ...,  cachename );
function value= get_cache_obj( obj, prop, varargin );
   global eidors_objects
   value= [];

   [objlist, cachename]= proc_obj_list( varargin{:} );

   try
      obj_id= obj.id;
   catch
      return; % Can't cache onto non-existant objects
   end

   if nargin==3
      for dep_obj = objlist{:}
         prop= [prop,'_', calc_obj_id(dep_obj)];
      end
   end

% if cachename is specified, then cache to that file, rather
%  than to the standard eidors_objects location
% TODO: fixthis - use ( ) for matlab > 6.0
   if ~isempty( cachename )
      test_for_cachdir;
      filename= [ eidors_objects.cachedir, '/' , prop, '_' cachename '.mat' ];
      if exist( filename, 'file')
         load(filename); %variable 'value' should be there
      end
   else
      try
   %     value= eidors_objects.( obj_id ).cache.( prop );
         value= eval(sprintf('eidors_objects.%s.cache.%s;',obj_id,prop));
      end
   end

function set_cache_obj( obj, prop, value, varargin )
   global eidors_objects
   [objlist, cachename]= proc_obj_list( varargin{:} );

   try
      obj_id= obj.id;
   catch
      return; % Can't cache onto non-existant objects
   end

   if nargin>=4 & isempty(cachename)
      for dep_obj = objlist{:}
         prop= [prop,'_', calc_obj_id(dep_obj)];
      end
   end

   if isempty(cachename)
   %  eidors_objects.( obj_id ).cache.( prop ) = value;
      eval(sprintf('eidors_objects.%s.cache.%s=value;', obj_id, prop));
   else
      filename= [ eidors_objects.cachedir, '/' , prop, '_' cachename '.mat' ];
      save(filename, 'value');
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
function obj_id= calc_obj_id( var )
   if isfield(var,'id')
      obj_id= var.id;
      return;
   end

   global eidors_objects;
   if ~isfield(eidors_objects,'hash_type')
       test_for_hashtypes;
   end

   if eidors_objects.hash_type==1; 
       % the obj_id does not depend on the cache, id or name
       try; var = rmfield(var,'cache'); end
       try; var = rmfield(var,'id');    end
       try; var = rmfield(var,'name');  end

       obj_id= eidors_var_id( var );
   elseif eidors_objects.hash_type > 1e6
%if hashing code is unavailable, then disable caching function
       obj_id= sprintf('id_%08d', eidors_objects.hash_type );
       eidors_objects.hash_type= eidors_objects.hash_type + 1;
   else
       error('hash_type value unrecognized');
   end

% Since we use a dynamically loaded function to test
% for hashtypes, all sorts of things can go wrong
%   - The compilers or libraries may not be available
%   - Compiler library versions may not match Matlab
%
% We test for existance of 'eidors_var_id' mex file
% 
% If nothing exists, then the best we can do is to
% disable hashing completely. We set hashtype to a
% number and increment it each time   
function test_for_hashtypes
   global eidors_objects; 
   if exist('eidors_var_id')==3 % MEX-file on MATLAB's search path
      eidors_objects.hash_type = 1;
   else
      eidors_objects.hash_type= 1e8+1; 
   end

% Test whether the cachedir field has been set. This is
%  where eidors will store cached calculations. If it has
%  not been set, then create it as 'eidors_cache' in the
%  current directory
function test_for_cachdir
   global eidors_objects; 
   if ~isfield(eidors_objects, 'cachedir')
      cachedir= 'eidors_cache';
      eidors_objects.cachedir= [pwd,'/',cachedir];

      if ~exist( eidors_objects.cachedir, 'dir')
% Now we need to ensure that cachedir exists. Because the
% STUPID!!! matlab has a completely useless and nonstandard
% mkdir function, we try this
         mkdir(pwd,cachedir); 
      end
   end

% Test whether a cachedir function has been provided
% (by testing whether the last entry is a string). If so,
% return it in objlist, otherwise return []
function [objlist, cachedir]= proc_obj_list( varargin );
   cachedir= [];
   if nargin==0
      objlist= {};
   elseif isstr(varargin{nargin})
      cachedir= varargin{nargin};
      objlist = varargin(1:nargin-1);
   else
      objlist = varargin(:);
   end
