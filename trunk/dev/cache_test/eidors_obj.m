function obj_id= eidors_obj(type,name, varargin )
% EIDORS_OBJ: maintains EIDORS internals
%
% USAGE: to get eidors_version
%     version = eidors_obj('eidors_version')
%
% USAGE: to get interpreter version:
%     version = eidors_obj('interpreter_version')
%
% USAGE: to get path to EIDORS:
%     path = eidors_obj('eidors_path');
%
% USAGE: to cache values (not recommended)
%          eidors_obj('set-cache',obj, cachename,value, [time])
%     obj= eidors_obj('get-cache',obj, cachename)
%
% this will get or set the values of cached properties of the object.
%
%    example: % set jacobian
%        eidors_obj('set-cache',cache_obj, 'jacobian', J);
%
%    example: % get jacobian or '[]' if not set
%        J= eidors_obj('get-cache',cache_obj, 'jacobian');
%
% It is recommended to combine in cache_obj the minimum set of variables on
% which the value to be cached depends.
%    example: % cache_obj for jacobian
%        cache_obj = {img.fwd_model.nodes, img.fwd_model.elems ...
%                     img.elem_data, img.fwd_model.jacobian}
%
% NOTE that rather than directly using eidors_obj to set and get cache, it 
% is recommended to use eidors_cache with a function_handle.
%
%   example:
%        J = eidors_cache(@calc_jacobian_adjoint,img,'jacobian');
%
% See also: EIDORS_CACHE

% (C) 2005-10 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin==0 || ~ischar(type)
   error('cannot call eidors_obj with no arguments');
end

switch type
   case 'set'
      obj_id= set_obj( name, varargin{:} );
   case 'get-cache'
      test_install
      obj_id = [];
      if status_check(varargin{1})
        obj_id= get_cache_obj( name, varargin{:} );
      end
      
   case 'set-cache'
      test_install
      obj_id= [];
      if status_check(varargin{1})
          set_cache_obj( name, varargin{:} );
      end

   case 'eidors_version'
      obj_id= '3.8+';  % Update for New eidors version
      
   case 'eidors_path'
      global eidors_objects
      obj_id = eidors_objects.eidors_path;

   case 'interpreter_version'
      obj_id= test_versions;
% TODO: Add these functions
%  case 'eidors_path'
%  case 'eidors_dev_path'
%  case 'eidors_cache_path'

   case 'cache_init'
      cache_init;
      
   otherwise
      test_install
      obj_id= new_obj( type, name, varargin{:} );
end

function ok = status_check(name)
ok = true;
switch cache_status
    case 0
        ok = false;
    case 0.5
        dbs = dbstack;
        if cache_status(dbs(3).name) == 0
            ok = false;
        end
        if strcmp(dbs(3).name,'cache_shorthand') && cache_status(dbs(5).name) == 0
            ok = false;
        end
        if cache_status(name) == 0
           ok = false;
        end
end

function on = debug_status_check
on = false;
switch debug_status
    case 1
        on = true;
    case 0.5
        dbs = dbstack;
        if debug_status(dbs(3).name) == 1
            on = true;
        end
end


function out = cache_status(fname)
    global eidors_objects;
    if nargin == 0
        try
            out = eidors_objects.cache_enable;
        catch
            out = 1;
        end
    else
        out = ~any(strcmp(eidors_objects.cache_disabled_on,fname));
    end

    
function out = debug_status(fname)   
   global eidors_objects;
   if nargin == 0
      try
         out = eidors_objects.debug_enable;
      catch
         out = 0;
      end
   else
      out = ~any(strcmp(eidors_objects.debug_enabled_on,fname));
   end
      
      
function test_install
  global eidors_objects;
  if isfield(eidors_objects,'max_cache_size'); return; end % OK
  error('EIDORS not correctly started. Did you do ">>run /path/to/eidors/startup"');

function verstr = test_versions;
      ver= version; ver(ver=='.')=' ';
      ver = sscanf(ver,'%f'); ver=ver(:);

      % Convert 7.8.1 to 7.008001
      verstr.ver = 1e-3.^(0:length(ver)-1) * ver(:); 

      if exist('OCTAVE_VERSION') == 5
         verstr.isoctave = 1;
      else
         verstr.isoctave = 0;
      end

      cptr = computer;
      if strcmp(cptr(end+(-1:0)), '64')
         verstr.is64bit = 1;
      else
         verstr.is64bit = 0;
      end

function obj = set_obj( obj, varargin );
   global eidors_objects;

%  eidors_objects.( obj_id ) = obj;
%  eidors_objects.( obj_id ).cache= []; %clear cache
%  

   for idx= 1:2:nargin-1
%     obj.( obj_id ).( varargin{idx} )= varargin{idx+1};
% for matlab 6.1 compatibility 
      eval(sprintf('obj.%s=varargin{%d};', varargin{idx},idx+1 ));
   end

% val= get_cache_obj( obj, prop, dep_obj1, dep_obj2, ...,  cachename );
function value= get_cache_obj( obj, prop )
   global eidors_objects
   value= [];

% We don't do this since cache directories aren't defined (yet)
%  [objlist, cachename]= proc_obj_list( varargin{:} );
   if ~isfield(eidors_objects, 'cache')
       cache_init;
       return
   end

   if isempty(eidors_objects.cache.meta), 
      return 
   end
   c = eidors_objects.cache.cols;
%    match = ismember(prop, eidors_objects.cache.meta(:,c.prop));
%    if any(match)
   obj_id= calc_obj_id( { obj, prop} ); % recalculate in case obj changed
      
% if cachename is specified, then cache to that file, rather
%  than to the standard eidors_objects location
% TODO: fixthis - use ( ) for matlab > 6.0
%  if ~isempty( cachename )
%     test_for_cachdir;
%     filename= [ eidors_objects.cachedir, '/' , prop, '_' cachename '.mat' ];
%     if exist( filename, 'file')
%        load(filename); %variable 'value' should be there
%     end
%  else
      try
         value= eidors_objects.cache.( obj_id );
         idx = find(strcmp(obj_id, eidors_objects.cache.meta(:,c.obj_id)), 1, 'first');
         eidors_objects.cache.meta{idx,c.time} = now;
         eidors_objects.cache.meta{idx,c.count} = eidors_objects.cache.meta{idx,c.count} + 1;
         eidors_objects.cache.meta{idx,c.score_eff} = calc_effort_score(...
            eidors_objects.cache.meta{idx,c.effort}, ...
            eidors_objects.cache.meta{idx,c.count}, ...
            eidors_objects.cache.meta{idx,c.prio});
%        value= eval(sprintf('eidors_objects.%s.cache.%s;',obj_id,prop));
%          check_size(obj_id, prop);
      catch err
         if ~strcmp(err.identifier,'MATLAB:nonExistentField')
            rethrow(err);
         end
      end
%  end
%    end

function set_cache_obj( obj, prop, value, time )
   global eidors_objects
   if ~cache_this( obj ) ; return ; end

   if nargin  < 4
      time = 1; % assume 1 sec
   end
   
   if ~isfield(eidors_objects,'cache');
      init_cache;
   end

   c = eidors_objects.cache.cols;
   
% Cache directories aren't defined, yet
%  [objlist, cachename]= proc_obj_list( varargin{:} );

   obj_id = calc_obj_id( {obj, prop} );
   
   prio = eidors_objects.cache_priority;

%  if isempty(cachename)

      ws = whos('value');
      row(c.obj_id)   = {obj_id};
      row(c.prop)     = {prop};
      row(c.size)     = {ws.bytes};
      row(c.score_sz) = {calc_size_score(ws.bytes)};
      row(c.effort)   = {time};
      row(c.prio)     = {prio};
      row(c.count)    = {1};
      row(c.score_eff)= {calc_effort_score(time, 1, prio)};
      row(c.time)     = {now};
      
      if isfield(eidors_objects.cache, obj_id)
         idx = find(strcmp(obj_id, eidors_objects.cache.meta(:,c.obj_id)));
         eidors_msg('@@ replaced cache object %s { %s }', obj_id, prop,4);
         eidors_objects.cache.size = eidors_objects.cache.size ...
                                    - eidors_objects.cache.meta{idx,c.size};
                                
      else
         idx = size(eidors_objects.cache.meta, 1) + 1;
      end
      eidors_objects.cache.meta(idx,:) = row;
      eidors_objects.cache.( obj_id ) = value;
      eidors_objects.cache.size = eidors_objects.cache.size + ws.bytes;
      check_size(obj_id, prop);
%  else
%     filename= [ eidors_objects.cachedir, '/' , prop, '_' cachename '.mat' ];
%     save(filename, 'value');
%  end

function cache_init
   global eidors_objects;
   
   eidors_objects.cache = struct;
   eidors_objects.cache.meta = cell(0);
   eidors_objects.cache.cols.obj_id    = 1;
   eidors_objects.cache.cols.prop      = 2;
   eidors_objects.cache.cols.time      = 3;
   eidors_objects.cache.cols.size      = 4;
   eidors_objects.cache.cols.score_sz  = 5;
   eidors_objects.cache.cols.effort    = 6;
   eidors_objects.cache.cols.prio      = 7;
   eidors_objects.cache.cols.count     = 8;
   eidors_objects.cache.cols.score_eff = 9;
   eidors_objects.cache.size           = 0;

function score = calc_size_score(sz)
   if iscell(sz)
      fn = @(b) round(10*log10(b / 1024));
      N = numel(sz);
      score = num2cell(cellfun(fn, sz));
   else
      score = round(10 * log10( sz / 1024));
   end
   
function score = calc_effort_score( time, counts, prios)
   if iscell(time)
       score = num2cell( round(10*log10( cell2mat(time) .* cell2mat(counts))) +...
                                cell2mat(prios));
   else
      score = round(10*log10(time .* counts)) + prios; 
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
   try 
      obj_id= eidors_var_id( var );
   catch
      global eidors_objects;
      if ~isfield(eidors_objects,'hash_type')
         eidors_objects.hash_type= 1e8+1; 
      end
%if hashing code is unavailable, then disable caching function
      obj_id= sprintf('id_%031d%08d', 0,eidors_objects.hash_type );
      eidors_objects.hash_type= eidors_objects.hash_type + 1;
   end

% Test whether the cachedir field has been set. This is
%  where eidors will store cached calculations. If it has
%  not been set, then create it as 'eidors_cache' in the
%  current directory
%
% NOTE: This code is not (yet) used
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
   elseif ischar(varargin{nargin})
      cachedir= varargin{nargin};
      objlist = varargin(1:nargin-1);
   else
      objlist = varargin(:);
   end

function retval= cache_this( obj )
% we choose not to cache data and images because this will
% tend to fill up the workspace.
% TODO: make the DONT_CACHE list use configuable
   if ~isstruct( obj); retval=1; return; end
   DONT_CACHE= {'data','image'};
   if any(strcmp( obj.type, DONT_CACHE));
      retval = 0;
   else
      retval = 1;
   end

function check_size( obj_id , prop )
   global eidors_objects;   
%    eidors_objects.( obj_id ).( prop ).last_used = now;

   max_memory= eidors_objects.max_cache_size;
%    ww= whos('eidors_objects');
   if eidors_objects.cache.size  > max_memory
      eidors_cache('clear_max',floor(max_memory*.75));
   end

