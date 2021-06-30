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
%   eidors_cache( 'list' )
%   eidors_cache( 'show_objs' )
%       - list all objects
%   eidors_cache(___, order)
%       - sort by specific order. Valid values are 
%         {time}, size, effort, count, rank
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
%        opt.cache_to_disk = true or '.'  
%        opt.cache_to_disk = 'path/to/cache/dir'
%            - Save the cached file to disk and reload it when called
%            - The file is named 'eidors_cache_id????.mat' or opt.fstr,'_id???'
%        opt.cache_on_ng_opt = true
%            - Cache also on the contents of the 'ng.opt' file in current working dir
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

if nargin==1 && ischar(command) && strcmp(command,'UNIT_TEST');
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


% In order to do caching, we need to have an existing
% function, as a string which exists as a function 
% handle.
if isa(command, 'function_handle') || ...
  (ischar(command) && ...
   any(exist(command) == [2,3])) % m or mex-file
   % Test if we're being asked to cache on anonymous function
   %   and issue error 
   if isa(command, 'function_handle')
      str = func2str(command);
      if str(1) == '@'
         error('Cannot cache anonymous functions');
      end
   end
   [varargout{1:nargout}] =  ...
           cache_shorthand(command, varargin{:});
   return
end

    
switch command
   case 'init'
      eidors_objects.cache_enable = 1;
      eidors_objects.cache_disabled_on = [];
      eidors_objects.cache_debug_enable = 0;
      eidors_objects.cache_debug_enabled_on = [];
      eidors_obj('cache_init');
      
   case 'clear_all'
      remove_objids

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
      if ischar(limit); limit= str2double(limit); end
         varargout{1} = varargout{1} + limit;
      end
      eidors_objects.cache_priority = varargout{1};

   case {'list', 'show_objs'} 
      if nargin == 2 
         cache_list(limit);
      else
         cache_list;
      end
      
      
   case 'clear_max'
      if ischar(limit); limit= str2double(limit); end
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

function cache_list (order)
   global eidors_objects;
   try
      meta = eidors_objects.cache.meta;
   catch
      return
   end
   if nargin == 0
       order = 'time';
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
   switch order
       case 'time'
           meta = mysortrows(meta,c.time); % sort by time
       case {'prop','name'}
           meta = mysortrows(meta,c.prop); % sort by name
       case 'rank'
           meta = mysortrows(meta,N+1); % sort by rank
       case 'size'
           meta = mysortrows(meta,-c.size); % sort by rank
       case 'effort'
           meta = mysortrows(meta,-c.effort);
       case 'count'
           meta = mysortrows(meta,-c.count);
       otherwise
           error('Unrecognized sort order');
   end
   
   meta = meta';
   fprintf('CACHE__:  Date+Time     bytes       Score_szPrio CountXEffort   Score_eff   #  obj_id ___\n');
   fprintf('%s b=%9.0d [%4d]  p=%02d t=%3dx%.2e [%4d] i=%4d: %s { %s }\n', ...
      meta{[c.time,c.size,c.score_sz,c.prio,c.count,c.effort,c.score_eff,N+1,c.obj_id, c.prop],:});

function priidx = get_cache_priority
   global eidors_objects;
   priidx = [];
   if isfield(eidors_objects.cache, 'meta') && isfield(eidors_objects.cache, 'cols')
      meta = eidors_objects.cache.meta;
      c = eidors_objects.cache.cols;
      [jnk, priidx] = mysortrows(meta,[-c.score_eff c.score_sz -c.time]);
      priidx(priidx) = 1:size(meta,1);
   end

      
   
function objid = clear_names_cache( name )
   objid=[];
   global eidors_objects;
   try
      c = eidors_objects.cache.cols;
      objid = find(strcmp(name, eidors_objects.cache.meta(:,c.prop)));
   end
   
   
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
      eidors_objects.cache.size = eidors_objects.cache.size - total_size;
   end
   
   eidors_msg('Removed %d objects with %d bytes from cache', ...
      N, total_size, 2 );
   
       
%   v1 = eidors_cache( @function_handle, {param1, param2, ...})
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
   cache_to_disk = false;
   if isfield(opt, 'cache_to_disk')
       cache_locn    = opt.cache_to_disk;
       if cache_locn
          cache_to_disk = true;
       end
       if cache_locn==true
          cache_locn = '.';
       end
       if isfield(opt,'fstr')
          cache_str = opt.fstr;
       else
          cache_str = 'eidors_cache';
       end
   end
   if ~isfield(opt, 'cache_on_ng_opt')
       opt.cache_on_ng_opt = false;
   end
   if opt.cache_on_ng_opt
     cache_obj{end+1} = cache_ng_opt_bytes();
   end

   [varargout,obj_id] = eidors_obj('get-cache', cache_obj, fstr );
   if length(varargout)==0 && cache_to_disk % if not in memory cache
       savename = [cache_locn,'/',cache_str,'_',obj_id,'.mat'];
       if exist(savename,'file')
          load(savename)
       end  
   end
   if numel(varargout) < nargout
      eidors_msg('@@ (Re)calculating %s',fstr, level_in);
      t0 = tic;
      clear varargout; % to ensure it still runs on older MLs (e.g. R2012b)
      [varargout{1:nargout}] =  ...
               feval(fhandle, args{:});
      t = toc(t0);
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',opt.boost_priority);
      end
      

      if cache_to_disk
         eidors_msg('@@ Caching to %s', savename, level_in+1);
         save(savename,'varargout','-V7'); % -V7 for octave support
      else
         eidors_obj('set-cache', cache_obj, fstr, varargout, t);
      end
      
      if isfield(opt,'boost_priority');
         eidors_cache('boost_priority',-opt.boost_priority);
      end
      return
   end
   eidors_msg('%s: Using cached value',fstr,level_out);

function bytes = cache_ng_opt_bytes
   if ~exist('ng.opt');
      bytes = []; return
   end
   fid = fopen('ng.opt','rb');
   bytes = uint8(fread(fid,[1,inf],'uint8'));
   fclose(fid);

% replace meshsizefilename with contents
   mszstr = 'options.meshsizefilename  ';
   ff = findstr(bytes,mszstr);
   if isempty(ff); return; end
   if length(ff)>1 
      eidors_msg('Unexpected bug. Check ng.opt writer',1);
   end
   bytescut(1) = ff(1)+length(mszstr);
   ff = findstr(bytes(bytescut(1):end),'.msz');
   if length(ff)==0;
      bytes = []; return % empty meshsizefilename. That's OK.
   end
   bytescut(2) = bytescut + ff(1) + 2;
   mszfile = bytes(bytescut(1):bytescut(2));
   fid = fopen(char(mszfile),'rb');
   if fid>0
       bytes2= uint8(fread(fid,[1,inf],'uint8'));
       fclose(fid);
   else
       bytes2 = [];
       eidors_msg('Warning. no MSZ file %s',mszfile,1);
   end
   bytes = [bytes(1:bytescut(1)-1), ...
            bytes2, ...
            bytes(bytescut(2)+3)];

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
   test_disk_caching
 
function [v1 v2] = test_function(a,b,c,d)
   v1 = rand(1);
   v2 = rand(1);
   
function [meta,idx] = mysortrows(meta, cols)
   if ~exist('OCTAVE_VERSION') % octave doesn't sort cells (valid jun 2019)
      [meta,idx] = sortrows(meta, cols);
   else
      metm = cell2mat(meta(:,3:end));
      cols = ( abs(cols) - 2 ) .* sign(cols);
      [~,idx] = sortrows(metm, cols);
      meta = meta(idx,:);
   end


function test_debug
   fprintf('\n\n************************\n        CACHE DEBUG: VERBOSE OUTPUT\n************************\n');
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
   fprintf('\n\n************************\n        CACHE DEBUG: DEBUG FINISHED\n************************\n');
   
function test_disk_caching
   tstdir = 'eidors_cache_test_dir1234321';
   [~]=rmdir(tstdir,'s');
   [~]=mkdir(tstdir);
   disp('============= DISK =============');
   opt = struct('cache_to_disk',tstdir);
   eidors_cache clear; global eidors_objects
   esz = length(fieldnames(eidors_objects.cache));
   unit_test_cmp('disk 00',length(dir([tstdir,'/*.mat'])),0)
   [v1]= eidors_cache( @sin, {1:100},opt);
   unit_test_cmp('disk 01',length(dir([tstdir,'/*.mat'])),1)
   unit_test_cmp('disk 02',0, ...
      length(fieldnames(eidors_objects.cache))-esz);
   obj_id = eidors_var_id({1:100,'sin'});
   cfname = [tstdir,'/eidors_cache_',obj_id,'.mat'];
   unit_test_cmp('disk 03',length(dir(cfname)),1); 
   loader = load(cfname);
   unit_test_cmp('disk 04',loader.varargout{1},v1);

   [v1]= eidors_cache( @sin, {1:100},opt);
   unit_test_cmp('disk 05',length(dir([tstdir,'/*.mat'])),1)
   unit_test_cmp('disk 06',0, ...
      length(fieldnames(eidors_objects.cache))-esz);
   loader = load(cfname);
   unit_test_cmp('disk 07',loader.varargout{1},v1);
   [~]=rmdir(tstdir,'s');
   disp('============= ~DISK =============');
   eidors_cache clear; global eidors_objects
   opt = struct('cache_to_disk',false);
   esz = length(fieldnames(eidors_objects.cache));
   [v1]= eidors_cache( @sin, {1:100},opt);
   unit_test_cmp('not disk 01',1, ...
      length(fieldnames(eidors_objects.cache))-esz);
   obj_id = eidors_var_id({1:100,'sin'});
   unit_test_cmp('not disk 02',1, ...
      length(fieldnames(eidors_objects.cache))-esz);
   unit_test_cmp('not disk 03',isfield(eidors_objects.cache,obj_id),1);
   unit_test_cmp('not disk 04',eidors_objects.cache.(obj_id),[v1]);

   [v1]= eidors_cache( @sin, {1:100},opt);
   unit_test_cmp('not disk 05',1, ...
      length(fieldnames(eidors_objects.cache))-esz);

   
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
