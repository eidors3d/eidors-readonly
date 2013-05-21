function varargout=eidors_cache( command, varargin )
% EIDORS_CACHE Control eidors caching
% Usage: eidors_cache( command, limit ) for cache management
% Usage: eidors_cache(@function_handle, {params1, ... }) for shorthand
%        chashing in m-files
%
% USAGE:
%   eidors_cache( 'clear_all' ) 
%   eidors_cache  clear
%
%   eidors_cache( 'clear_old'. timestamp );
%   eidors_cache( 'clear_new', timestamp );
%      - clear all variables older (or newer) than timestamp
%      - example:
%         time_now= now;
%         lots_of_eidors_calcs % don't need to cache this stuff
%         eidors_cache('clear_new',time_now);
%
%   eidors_cache( 'clear_max', memory_in_bytes );
%      - clear cache so it has less than memory_in_bytes
%
%   eidors_cache( 'cache_size', memory_in_bytes );
%      - set max cache size to be memory_in_bytes
%      - without 2nd arg will return current cache_size
%
%   eidors_cache( 'clear_name', cache_name )
%      - eg. eidors_cache( 'clear_name', 'inv_solve_diff_GN_one_step')
%      - clear all variables with name 
%
%   eidors_cache( 'clear_model_library' );
%      - clear the eidors model library directory
%      - NOTE: this function is experimental. 
%      - TODO: add a way to specify what to clear
%  
%   eidors_cache( 'boost_priority', value)
%      - modify the priority of the next cached items
%      - low priority variables will be deleted first when memory is tight
%
%   eidors_cache( 'off' )
%   eidors_cache( 'disable' )
%   eidors_cache( 'off', 'function_name' )
%      - new values are not added to cache
%      - requests for cached values return []
%      - if specified, only applies to a specific function
%
%   eidors_cache( 'on' )
%   eidors_cache( 'enable' )
%   eidors_cache( 'on', 'function_name' )
%      - re-enables caching
%      - if specified, only applies to a specific function
%
%   eidors_cache( 'status' )
%      - queries the caching status
%      - 0  : all off
%      - 1  : all on
%      - 0.5: off for some functions
%
%   eidors_cache( 'debug_on' )
%   eidors_cache( 'debug_on', 'function_name' )
%      - enables debug output
%      - if specified, only applies to a specific function
%      - will print a message everytime an object is removed from cache
%
%   eidors_cache( 'debug_off' )
%   eidors_cache( 'debug_off', 'function_name' );
%      - disables debug output
%      - if specified, only applies to a specific function
%
%   eidors_cache( 'debug_status' )
%      - queries debug status
%      - output analogous to cache status above
%
%   v1 = eidors_cache( @function_handle, {param1, param2, ...})
%   [v1, v2, ...] = eidors_cache( @function_handle, {param1, param2, ...})
%   [v1, v2, ...] = eidors_cache( ... , opt)
%      Shorthand interface for caching the output of a specific function.
%      Specify all arguments to the function as a cell array. They will all
%      be used for as a caching object by default. The following options
%      are available:
%        opt.cache_obj
%            - a single value, or a cell array of values to cache on, 
%              rather than all the inputs (e.g. not all values of an input
%              struct may be used in the function)
%        opt.boost_priority
%            - priority boost to use for that function
%
%   eidors_cache( 'cache_path' )
%   eidors_cache( 'cache_path', '/path/to/cache/path' )
%       - get and set cache_path, a path to a writable
%           directory in which eidors can store files
%
% See also EIDORS_OBJ

% (C) 2005-2013 Andy Adler and Bartlomiej Grychtol.
% License: GPL version 2
% $Id$

% Comments
% Want to clear specific structures
%      to clear old variables
%      to clear specific parts of structures

if nargin==1 && isstr(command) && strcmp(command,'UNIT_TEST');
      do_unit_test; return; end

global eidors_objects;
if nargin<1
   fprintf('EIDORS_CACHE: current max memory = %.0f MB\n', ...
         eidors_objects.max_cache_size/(1024*1024)); 
   ww= whos('eidors_objects');
   fprintf('EIDORS_CACHE: cache memory used = %.0f MB\n', ...
         ww.bytes/(1024*1024)); 
   fprintf('EIDORS_CACHE: current priority = %d\n', ...
         eidors_objects.cache_priority); 
   return;
elseif nargin > 1 
   limit = varargin{1};
end

if isa(command, 'function_handle')
   output = mk_varargout_str(nargout);
   eval(sprintf('%s = %s', output, 'cache_shorthand(command, varargin{:});'));
   return
end

    
switch command
   case 'init'
      eidors_objects.cache_enable = 1;
      eidors_objects.cache_disabled_on = [];
      eidors_objects.debug_enable = 0;
      eidors_objects.debug_enabled_on = [];
      
   case {'clear_all','clear'}
      [objid, times, sizes, prios, priidx] = get_names_times;
      remove_objids( objid, sizes,  1:length(sizes) );

   case 'cache_size'
      if nargin==2
      if ischar(limit); limit= str2num(limit); end
         eidors_objects.max_cache_size = limit;
      else
         varargout{1}= eidors_objects.max_cache_size;
      end

   case 'cache_path'
      if nargin == 1
         varargout{1}= eidors_objects.cache_path;
      else
         eidors_objects.cache_path = varargin{1};
      end
      
   case {'disable' 'off'}
       if nargin == 1
           eidors_objects.cache_enable = 0;
           eidors_objects.cache_disabled_on = {};
       else
           eidors_objects.cache_enable = 0.5;
           if isfield(eidors_objects,'cache_disabled_on')
            if ~any(strcmp(eidors_objects.cache_disabled_on, limit))
                eidors_objects.cache_disabled_on = [...
                    eidors_objects.cache_disabled_on; {limit}];
            end
           else
               eidors_objects.cache_disabled_on =  {limit};
           end
       end
   case {'enable' 'on'}
       if nargin == 1
           eidors_objects.cache_enable = 1;
           eidors_objects.cache_disabled_on = {};
       else
           if isfield(eidors_objects,'cache_disabled_on')
               idx = strcmp(eidors_objects.cache_disabled_on, limit);
               eidors_objects.cache_disabled_on(idx) = [];
           else 
               eidors_objects.cache_disabled_on = [];
           end
           if isempty(eidors_objects.cache_disabled_on)
               eidors_objects.cache_enable = 1;
           end
       end
   case 'status'
      if nargin == 1
         try
            varargout{1} = eidors_objects.cache_enable;
         catch
            varargout{1} = 1;
         end
      else
         if isfield(eidors_objects,'cache_disabled_on')
            idx = strcmp(eidors_objects.cache_disabled_on, limit);
            varargout{1} = double(~any(idx));
         end
      end
   case 'debug_status'
      if nargin == 1
         varargout{1} = eidors_objects.debug_enable;
      else
         if isfield(eidors_objects,'debug_enabled_on')
            idx = strcmp(eidors_objects.debug_enabled_on, limit);
            varargout{1} = double(any(idx)) | eidors_objects.debug_enable==1;
         end
      end
   case 'debug_on'
       if nargin == 1
           eidors_objects.debug_enable = 1;
           eidors_objects.debug_enabled_on = {};
       else
           eidors_objects.debug_enable = 0.5;
           if isfield(eidors_objects,'debug_enabled_on')
            if ~any(strcmp(eidors_objects.debug_enabled_on, limit))
                eidors_objects.debug_enabled_on = [...
                    eidors_objects.debug_enabled_on; {limit}];
            end
           else
               eidors_objects.debug_enabled_on =  {limit};
           end
       end
   case 'debug_off'
      if nargin == 1
         eidors_objects.debug_enable = 0;
         eidors_objects.debug_enabled_on = {};
      else
         if isfield(eidors_objects,'debug_enabled_on')
            idx = strcmp(eidors_objects.debug_enabled_on, limit);
            eidors_objects.debug_enabled_on(idx) = [];
         else
            eidors_objects.debug_enabled_on = [];
         end
         if isempty(eidors_objects.debug_enabled_on)
            eidors_objects.debug_enable = 0;
         end
      end
   case 'boost_priority'
      try
         varargout{1}= eidors_objects.cache_priority;
      catch
         varargout{1}= 0; % default priority
      end
      if nargin==2
      if ischar(limit); limit= str2num(limit); end
         varargout{1} = varargout{1} + limit;
      end
      eidors_objects.cache_priority = varargout{1};

   case 'show_objs'
      [objid, times, sizes, prios, priidx, names] = get_names_times;
      for i=1:length(times)
         fprintf('t=%9.7f b=%9.0d p=%02d, i=%03d: %s { ', rem(times(i),1), ... %today
             sizes(i), prios(i), priidx(i), objid{i} );
         fprintf('%s ', names{i}{:});
         fprintf('}\n');
      end

   case 'clear_max'
      if ischar(limit); limit= str2num(limit); end
% This will remove just in order of priidx.
%     remove_objids( objid, sizes,  find(cumsum(sizes) > limit) );
% Remove in order of time + priority
      [objid, times, sizes, prios, priidx] = get_names_times;
      tot=     cumsum(sizes(priidx)); 
      remove = find(tot > limit);
      rmidx=   priidx(remove);
      remove_objids( objid, sizes,  rmidx);

   case 'clear_old'
      if ischar(limit); limit= str2num(limit); end
      [objid, times, sizes, prios, priidx] = get_names_times;
      remove_objids( objid, sizes,  ...
        find(times < limit) );

   case 'clear_new'
      if ischar(limit); limit= str2num(limit); end
      [objid, times, sizes, prios, priidx] = get_names_times;
      remove_objids( objid, sizes,  ...
        find(times > limit) );

   case 'clear_model_library'
      %TODO: add ways to select what to delete
      delete([eidors_objects.model_cache,'/*.mat']);

   case 'clear_name'
      [objid, sizes, clearidx] = clear_names_cache( limit );
      remove_objids( objid, sizes, clearidx);
   
   case 'dump'
      varargout{1} = eidors_objects;
      
   case 'load'
      eidors_objects = limit;
      
   otherwise
      error('command %s not understood',command);
end

function [objid, sizes, clearidx] = clear_names_cache( name );
   objid={}; times=[]; clearidx= [];
   global eidors_objects;
   if isempty(eidors_objects); return; end
   idx=1;
   for fn= fieldnames(eidors_objects)'
      fn1= fn{1};
      if fn1(1:3)== 'id_'
         objid{idx}= fn1;
         obj = getfield(eidors_objects, fn1 );
         ww= whos('obj');
         sizes(idx) = ww.bytes;

         cache= getfield(obj,'cache');
         if isfield(getfield(obj,'cache'), name);
           clearidx(end+1) = idx;
         end

         idx= idx+1;
      end
   end

% priidx is priority index, where 1 prio is 0.1 of a day
function [objid, times, sizes, prios, priidx, names] = get_names_times;
   objid={}; times=[]; sizes=[]; pri=[]; prios=[]; priidx=[];names={};
   global eidors_objects;
   if isempty(eidors_objects); return; end
   idx=1;
   for fn= fieldnames(eidors_objects)'
      fn1= fn{1};
      if fn1(1:3)== 'id_'
         objid{idx}= fn1;

         obj = getfield(eidors_objects, fn1 );
         try
            times(idx) = obj.last_used;
         catch
            times(idx) = 0; % old
         end
         try
            prios(idx) = obj.priority;
         catch
            prios(idx) = 0; % default
         end
         names(idx) = {fieldnames(obj.cache)};
         ww= whos('obj');
         sizes(idx) = ww.bytes;

         idx= idx+1;
      end
   end

   % sort from recent (high) to old (low)
   [times, idx] = sort(-times);
   times= -times;
   objid= objid(idx);
   sizes= sizes(idx);
   prios= prios(idx);
   names= names(idx);
   % sort from high to low priority and recent (high) to old (low)
   [jnk, priidx] = sortrows([-prios(:), -times(:)]);

function remove_objids( objid, sizes, idx)
   global eidors_objects;
   
   switch eidors_cache('debug_status')
      case 1
         for i = 1:length(idx)
            debug_msg(objid{idx(i)},'removed');
         end
      case 0.5
         for i = 1:length(idx)
            flds = fieldnames(eidors_objects.(objid{idx(i)}).cache);
            for j = 1:numel(flds)
               if eidors_cache('debug_status',flds)
                  debug_msg(objid{idx(i)}, 'removed');
                  break
               end
            end
         end
   end
   
   eidors_objects= rmfield( eidors_objects, objid( idx ) );
   eidors_msg('eidors_cache: removing %d objects with %d bytes', ...
          length(idx), sum(sizes(idx)), 3 );

       
function varargout = cache_shorthand(fhandle, varargin)
% Will cache on all function inputs, unless opt.cache_obj is specified
   if nargin >2
      opt = varargin{2};
   else
      opt = struct;
   end
   if isfield(opt, 'cache_obj');
      cache_obj = opt.cache_obj;
      if ~iscell(cache_obj)
         cache_obj = {cache_obj};
      end
   else
      cache_obj = varargin{1};
   end
   
   fstr = func2str(fhandle);
   
   varargout = eidors_obj('get-cache', cache_obj, fstr );
   if numel(varargout) < nargout
      eidors_msg('@@@ (Re)calculating %s',fstr, 1);
      output = mk_varargout_str(nargout);
      varargout = cell(0);
      eval(sprintf('%s = %s', output, 'feval(fhandle,cache_obj{:});'));
      
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',opt.boost_priority);
      end
      
      eidors_obj('set-cache', cache_obj, fstr, varargout);
      
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',-opt.boost_priority);
      end
      return
   end
   eidors_msg('%s: Using cached value',fstr,2);

function output = mk_varargout_str(N)
output = '[';
for i = 1:N
   output = [ output sprintf('varargout{%d} ',i)];
end
output = [ output ']' ];

function debug_msg(id,action)
global eidors_objects;
name = fieldnames(eidors_objects.(id).cache);
fprintf('EIDORS_CACHE: %s %s {', action, id);
fprintf(' %s', name{:});
fprintf(' }\n');
% dbstack could be useful too

function do_unit_test

   ll= eidors_msg('log_level');
   eidors_msg('log_level',5);
   eidors_cache
   eidors_cache('clear_all');
   eidors_cache
   eidors_obj('set-cache', rand(1) , 't1', rand(2e3));
   eidors_obj('set-cache', rand(1) , 't2', rand(2e3));
   eidors_obj('set-cache', rand(1) , 't3', rand(2e3));
   eidors_cache('clear_name','t3');
   eidors_cache
   eidors_cache('clear_max', 34e6);
   eidors_cache
   eidors_cache('boost_priority', 1);
   eidors_cache
   eidors_cache('boost_priority', -1);
   eidors_cache
   [v1] = eidors_cache(@test_function,{3,4});
   [v2] = eidors_cache(@test_function,{3,4});
   unit_test_cmp('shorthand 1 param:',v1,v2);
   [v3 v4] = eidors_cache(@test_function,{3,4});
   unit_test_cmp('expect fail', v3, v4);
   [v5 v6 ] = eidors_cache(@test_function,{3,4});
   unit_test_cmp('shorthand 2 params:',v4, v6);
   [v5 v6 ] = eidors_cache(@test_function,{3,4, 5}); %this should re-calc
   opt.cache_obj = 5;
   [v5 v6 ] = eidors_cache(@test_function,{3,4, 5}, opt); %this should re-calc
   [v7 v8 ] = eidors_cache(@test_function,{1,2, 3}, opt); %this should NOT
   unit_test_cmp('shorthand cache_obj:',v6, v8);
   eidors_cache clear_all
   opt = struct;
   [v7 v8 ] = eidors_cache(@test_function,{1,2, 3}, opt); %this should NOT
   opt.boost_priority = 2;
   [v7 v8 ] = eidors_cache(@test_function,{3,4, 5}, opt); %this should NOT
   eidors_cache show_objs
   eidors_cache
   eidors_msg('log_level',ll);
   
   test_debug
   
function [v1 v2] = test_function(a,b,c,d)
   v1 = rand(1);
   v2 = rand(1);
   
function test_debug
   eidors_cache clear
   eidors_obj('set-cache',{5}, 'test1',50);
   eidors_obj('set-cache',{5}, 'test2',500);
   eidors_obj('set-cache',{10}, 'test1',100);
   eidors_cache show_objs
   eidors_cache debug_off
   eidors_cache debug_on test2
   eidors_cache clear_name test2
   eidors_cache show_objs
%    eidors_obj('set-cache',{5}, 'test1',50);
   eidors_obj('set-cache',{5}, 'test2',500);
   eidors_cache debug_off
   eidors_cache debug_on test1
   eidors_cache clear_name test2
   eidors_cache('clear_max',0)
   eidors_cache('show_objs')
   eidors_cache debug_off
   
   
   
%    obj= eidors_obj('get-cache',{5}, 'test1');
