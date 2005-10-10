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
%
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

% $Id: eidors_obj.m,v 1.24 2005-10-10 03:12:55 aadler Exp $
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
      
   obj.id= obj_id;
   if any(strcmp( obj.type, DONT_CACHE)); return; end 

% set the obj_id into the obj after calculating the hash value
   if ~isfield(eidors_objects,obj_id)
      eval(sprintf('eidors_objects.%s=obj;',obj_id));
   end

function val= calc_or_cache( obj, funcname, varargin )

   val = get_cache_obj( obj, funcname, varargin{:} );
   if ~isempty( val );
      eidors_msg([funcname, ': using cached value'], 2);
   else
      [objlist, cachename]= proc_obj_list( varargin{:} );
      val = feval( funcname, obj, objlist{:} );
      set_cache_obj( obj, funcname, val, varargin{:} );
      eidors_msg([funcname, ': setting cached value'], 2);
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

   % loop through extra args and fill to prop or cachename
   if nargin>=3 & isempty(cachename)
      for dep_obj = objlist{:}
          try
              prop= [prop,'_', dep_obj.id];
          catch
              prop= [prop,'_', calc_obj_id(dep_obj)];
          end
          
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
          try
              prop= [prop,'_', dep_obj.id];
          catch
              prop= [prop,'_', calc_obj_id(dep_obj)];
          end
      end
   end

   if isempty(cachename)
   %  eidors_objects.( obj_id ).cache.( prop ) = value;
      eval(sprintf('eidors_objects.%s.cache.%s=value;', obj_id, prop));
   else
      filename= [ eidors_objects.cachedir, '/' , prop, '_' cachename '.mat' ];
      save(filename, 'value');
keyboard
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
   if ~isfield(eidors_objects,'hash_type')
       test_for_hashtypes;
   end

   if eidors_objects.hash_type==1; 
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
   elseif eidors_objects.hash_type > 1e6
%if hashing code is unavailable, then disable caching function
       if isfield(var,'id')
          obj_id= var.id;
       else
           obj_id= sprintf('id_%08d', eidors_objects.hash_type );
           eidors_objects.hash_type= eidors_objects.hash_type + 1;
       end
   else
       error('hash_type value unrecognized');
   end
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
      eidors_objects.hash_type= 1e8+1; 
   end

% Test whether the cachedir field has been set. This is
%  where eidors will store cached calculations. If it has
%  not been set, then create it as 'eidors_cache' in the
%  current directory
function test_for_cachdir
   global eidors_objects; 
   if ~isfield(eidors_objects, 'cachedir')
      eidors_objects.cachedir= [pwd,'/eidors_cache'];

% Now we need to ensure that cachedir exists. Because the
% STUPID!!! matlab has a completely useless and nonstandard
% mkdir function, we shell out to the OS to do this for us.
% Luckly, most OSes provide a fairly reasonable mkdir function
      system(['mkdir ''',eidors_objects.cachedir,'']);
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
