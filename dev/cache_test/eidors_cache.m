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
%        opt.fstr
%            - name to use for the cached result
%        opt.log_level
%            - message level to use with eidors_msg. By default, it's 4 for
%              'setting cache' and 3 for 'using cached value'
%
%   eidors_cache( 'cache_path' )
%   eidors_cache( 'cache_path', '/path/to/cache/path' )
%       - get and set cache_path, a path to a writable
%           directory in which eidors can store files
%
%   eidors_cache( 'eidors_path' )
%   eidors_cache( 'eidors_path', '/path/to/eidors/' )
%       - /path/to/eidors is the path in which eidors_startup.m is found
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
      str = func2str(command);
      if str(1) == '@'
         error('Cannot cache anonymous functions');
      end
end

if isa(command, 'function_handle') || (ischar(command) && exist(command) == 2) % an m-file
   output = mk_varargout_str(nargout);
   eval(sprintf('%s = %s', output, 'cache_shorthand(command, varargin{:});'));
   return
end

    
switch command
   case 'init'
      eidors_objects.cache_enable = 1;
      eidors_objects.cache_disabled_on = [];
      eidors_objects.cache_debug_enable = 0;
      eidors_objects.cache_debug_enabled_on = [];
      
   case 'clear_all'
      remove_objids
%       [objid, times, sizes, prios, priidx, names] = get_names_times;
%       remove_objids( objid, names, sizes);

   case 'clear'
     switch nargin
       case 2
        eidors_cache('clear_name',limit);
       case 1
        eidors_cache('clear_all');
       otherwise
         error('Wrong number of inputs');
     end
      
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

   case 'eidors_path'
      if nargin == 1
         varargout{1}= eidors_objects.eidors_path;
      else
         eidors_objects.eidors_path = varargin{1};
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
         varargout{1} = eidors_objects.cache_debug_enable;
      else
         if isfield(eidors_objects,'cache_debug_enabled_on')
            idx = ismember(limit,eidors_objects.cache_debug_enabled_on);
            varargout{1} = idx | eidors_objects.cache_debug_enable==1;
         end
      end
   case 'debug_on'
       if nargin == 1
           eidors_objects.cache_debug_enable = 1;
           eidors_objects.cache_debug_enabled_on = {};
       else
           eidors_objects.cache_debug_enable = 0.5;
           if isfield(eidors_objects,'cache_debug_enabled_on')
            if ~any(strcmp(eidors_objects.cache_debug_enabled_on, limit))
                eidors_objects.cache_debug_enabled_on = [...
                    eidors_objects.cache_debug_enabled_on; {limit}];
            end
           else
               eidors_objects.cache_debug_enabled_on =  {limit};
           end
       end
   case 'debug_off'
      if nargin == 1
         eidors_objects.cache_debug_enable = 0;
         eidors_objects.cache_debug_enabled_on = {};
      else
         if isfield(eidors_objects,'cache_debug_enabled_on')
            idx = strcmp(eidors_objects.cache_debug_enabled_on, limit);
            eidors_objects.cache_debug_enabled_on(idx) = [];
         else
            eidors_objects.cache_debug_enabled_on = [];
         end
         if isempty(eidors_objects.cache_debug_enabled_on)
            eidors_objects.cache_debug_enable = 0;
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

   case {'list', 'show_objs'}
      cache_list
      
   case 'clear_max'
      if ischar(limit); limit= str2num(limit); end
      try
         c = eidors_objects.cache.cols;
      catch 
         return
      end
      
      priidx = get_cache_priority;
      [jnk, idx] = sort(priidx);
      tot=     cumsum([eidors_objects.cache.meta{idx,c.size}]); 
      rmidx  = idx(tot > limit);
      remove_objids( rmidx );

   case 'clear_old'
      if ischar(limit); limit= str2num(limit); end
      try
         c = eidors_objects.cache.cols;
      catch
         return
      end
      idx = find([eidors_objects.cache.meta{:,c.time}] < limit);
      remove_objids( idx );

   case 'clear_new'
      if ischar(limit); limit= str2num(limit); end
      try
         c = eidors_objects.cache.cols;
      catch
         return
      end
      idx = find([eidors_objects.cache.meta{:,c.time}] > limit);
      remove_objids( idx );

   case 'clear_model_library'
      %TODO: add ways to select what to delete
      delete([eidors_objects.model_cache,'/*.mat']);

   case 'clear_name'
      idx = clear_names_cache( limit );
      remove_objids( idx );
      
   case 'dump'
      varargout{1} = eidors_objects;
      
   case 'load'
      eidors_objects = limit;
      
   otherwise
      error('command %s not understood',command);
end

function cache_list
   global eidors_objects;
   try
      meta = eidors_objects.cache.meta;
   catch
      return
   end
   c = eidors_objects.cache.cols;
   if isempty(meta)
      fprintf('No objects in cache\n');
      return
   end
   meta(:,c.time) = cellstr(datestr([meta{:,c.time}],'yyyy-mm-dd HH:MM:SS.FFF'));
   N = size(meta,2);
   %    [jnk, priidx] = sortrows(meta(:,[c.score_eff c.score_sz c.time]),[-1 2 -3]);
   
   meta(:,N+1) = num2cell(get_cache_priority);
   
   meta = sortrows(meta,c.time); % sort by time
   
   meta = meta';
   fprintf('%s b=%9.0d [%4d]  p=%02d t=%3dx%.2e [%4d] i=%4d: %s { %s }\n', ...
      meta{[c.time,c.size,c.score_sz,c.prio,c.count,c.effort,c.score_eff,N+1,c.obj_id, c.prop],:});

function priidx = get_cache_priority
   global eidors_objects;
   priidx = [];
   try
      meta = eidors_objects.cache.meta;
      c = eidors_objects.cache.cols;
   end
   [jnk, priidx] = sortrows(meta,[-c.score_eff c.score_sz -c.time]);
   priidx(priidx) = 1:size(meta,1);
      
   
function objid = clear_names_cache( name )
   objid=[];
   global eidors_objects;
   try
      c = eidors_objects.cache.cols;
      objid = find(strcmp(name, eidors_objects.cache.meta(:,c.prop)));
   end
   
% priidx is priority index, where 1 prio is 0.1 of a day
function [objid, times, sizes, prios, priidx, names, effort] = get_names_times;
   objid={}; times=[]; sizes=[]; pri=[]; prios=[]; priidx=[];names={};
   effort=[];
   global eidors_objects;
   if isempty(eidors_objects); return; end
   idx=1;
   for fn= fieldnames(eidors_objects)'
      fn1= fn{1};
      if fn1(1:3)== 'id_'
         nms = fieldnames(eidors_objects.(fn1));
         for n = 1:numel(nms);
            objid{idx} = fn1;
            names{idx} = nms{n};
            times(idx) = 0; prios(idx) = 0; 
            effort(idx) = 10;
            obj = getfield(eidors_objects.(fn1), nms{n} );
            try times(idx) = obj.last_used; end
            try prios(idx) = obj.priority; end
            effort(idx) = obj.time_spent;
            ww= whos('obj');
            sizes(idx) = ww.bytes;
            
            idx= idx+1;
         end
      end
   end

   % sort from recent (high) to old (low)
   [times, idx] = sort(-times);
   times= -times;
   objid= objid(idx);
   sizes= sizes(idx);
   prios= prios(idx);
   names= names(idx);
   effort= effort(idx);
   
   priidx = calc_cache_priority(prios, times, effort, sizes);

   
function priidx = calc_cache_priority(prios, times, effort, sizes)
   % sort from high to low priority and recent (high) to old (low)
   % old way
%    [jnk, priidx] = sortrows([-prios(:), -times(:)]);
   % new way 
   effort_score = 2.^(round(10*log10(effort)) + prios) .* 2.^prios;
   size_score   = round(10*log10(sizes / 1024));
   [jnk, priidx] = sortrows([-effort_score(:), size_score(:), -times(:)]);
   priidx(priidx) = 1:length(sizes);

   
function remove_objids(idx, names, sizes)
   global eidors_objects;
   try
   c = eidors_objects.cache.cols;
   catch
      eidors_obj('cache_init'); 
      return % nothing else to do
   end
   if nargin == 0
      idx = 1:size(eidors_objects.cache.meta,1);
   end
   if isempty(idx)
      return
   end
   switch eidors_cache('debug_status')
      case 1
         debug_msg(  eidors_objects.cache.meta(idx,c.obj_id), ...
                     eidors_objects.cache.meta(idx,c.prop), 'removed');
      case 0.5
         db = eidors_cache('debug_status', eidors_objects.cache.meta(idx,c.prop));
         debug_msg(  eidors_objects.cache.meta(idx(db),c.obj_id), ...
                     eidors_objects.cache.meta(idx(db),c.prop), 'removed');
   end
   total_size = sum(cell2mat(eidors_objects.cache.meta(idx,c.size)));   
   N = numel(idx);
   if numel(idx) == size(eidors_objects.cache.meta,1)
      eidors_objects = rmfield(eidors_objects,'cache');
      eidors_obj('cache_init');
   else
      eidors_objects.cache = rmfield(eidors_objects.cache, ...
         eidors_objects.cache.meta(idx,c.obj_id));
      eidors_objects.cache.meta(idx,:) = [];
   end
   
   eidors_msg('Removed %d objects with %d bytes from cache', ...
      N, total_size, 2 );
   
function remove_objids_old( objid, names, sizes)
   global eidors_objects;
   N = length(objid);
   if ~iscell(names) 
      names = {names};
   end
   if numel(names) == 1
      names(1:N) = names;
   end
   switch eidors_cache('debug_status')
      case 1
         for i = 1:N
            debug_msg(objid{i},names{i},'removed');
         end
      case 0.5
         for i = 1:N
            if eidors_cache('debug_status',names{i})
               debug_msg(objid{i}, names{i}, 'removed');
               break
            end
         end
   end
   if nargin < 3
      sizes = zeros(1,N);
      for i = 1:N
         obj = getfield(eidors_objects.(objid{i}), names{i});
         ww = whos('obj');
         sizes(i) = ww.bytes;
      end
   end
   for i = 1:N
      eidors_objects.(objid{i})= rmfield( eidors_objects.(objid{i}), names{i} );
      if isempty(fieldnames(eidors_objects.(objid{i})))
         eidors_objects = rmfield( eidors_objects, objid{i});
      end
   end
   eidors_msg('Removed %d objects with %d bytes from cache', ...
          N, sum(sizes), 2 );

       
function varargout = cache_shorthand(fhandle, varargin)
% Will cache on all function inputs, unless opt.cache_obj is specified
   args = varargin{1};
   if ~iscell(args)
      args = {args};
   end
   
   if nargin >2
      opt = varargin{2};
   else
      opt = struct;
   end
   if ischar(opt)
      fstr = opt; clear opt;
      opt.fstr = fstr;
   end
   if isfield(opt, 'cache_obj');
      cache_obj = opt.cache_obj;
      if ~iscell(cache_obj)
         cache_obj = {cache_obj};
      end
   else
      cache_obj = args;
   end
   try
      fstr = opt.fstr;
   catch
      fstr = func2str(fhandle);
   end
   if isfield(opt, 'log_level')
       level_in = opt.log_level;
       level_out = opt.log_level;
   else
       level_in = 4;
       level_out = 3;
   end
   
   varargout = eidors_obj('get-cache', cache_obj, fstr );
   if numel(varargout) < nargout
      eidors_msg('@@ (Re)calculating %s',fstr, level_in);
      output = mk_varargout_str(nargout);
      varargout = cell(0);
      t0 = tic;
      eval(sprintf('%s = %s', output, 'feval(fhandle,args{:});'));
      t = toc(t0);
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',opt.boost_priority);
      end
      
      eidors_obj('set-cache', cache_obj, fstr, varargout, t);
      
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',-opt.boost_priority);
      end
      return
   end
   eidors_msg('%s: Using cached value',fstr,level_out);

function output = mk_varargout_str(N)
output = '[';
for i = 1:N
   output = [ output sprintf('varargout{%d} ',i)];
end
output = [ output ']' ];

function debug_msg(id,name,action)
global eidors_objects;
if nargin < 3 
   action = name;
   name = fieldnames(eidors_objects.(id));
end
if isempty(id), return, end
if ~iscell(name) name = {name}; end
if ~iscell(id) id = {id}; end

str = sprintf('EIDORS_CACHE: %s %%s { %%s }\\n', action);
arr = [id, name]'; 

fprintf(str, arr{:});
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
   eidors_cache
   eidors_cache list
   eidors_cache('clear_name','t3');
   eidors_cache
   eidors_cache list
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
   unit_test_cmp('Expect Fail', v3, v4,-inf);
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
   try
      eidors_cache(@(x) x^2, 3);
   catch
      eidors_msg('Error on anonymous function: correct',2);
   end
   test_debug
   test_priority
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
   
   
function test_priority
   eidors_cache clear
   eidors_obj('set-cache',{1}, 'slow_small', zeros(10) ,10); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('set-cache',{1}, 'slow_new',   zeros(100),10); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('set-cache',{1}, 'slow_big1',  zeros(100),10); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('set-cache',{1}, 'slow_big2',  zeros(100),20); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('set-cache',{1}, 'fast_small', zeros(2)  , 1); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('set-cache',{1}, 'fast_big',   zeros(100), 1); pause(.1)
%    eidors_cache list
%    fprintf('\n');
   eidors_obj('get-cache',{1}, 'slow_new'); 
   
   eidors_cache list

   
%    obj= eidors_obj('get-cache',{5}, 'test1');
