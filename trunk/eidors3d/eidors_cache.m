function retval=eidors_cache( command, limit )
% Control eidors_caching
% Usage: eidors_cache( command, limit )
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
%   eidors_cache( 'boost_priority', value)
%      - modify the priority of the next cached items
%      - higher priority variables will be deleted less when memory is tight  

% (C) 2005 Andy Adler. Licensed under GPL version 2
% $Id: eidors_cache.m,v 1.17 2007-08-30 03:37:04 aadler Exp $

% Comments
% Want to clear specific structures
%      to clear old variables
%      to clear specific parts of structures

if nargin<1
    error('eidors_cache requires at least one argument');
elseif nargin>=2
   if isstr(limit); limit= str2num(limit); end
end

global eidors_objects;

[objid, times, sizes, prios, priidx] = get_names_times;

switch command
   case {'clear_all','clear'}
      remove_objids( objid, sizes,  1:length(sizes) );

   case 'cache_size'
      if nargin==2
         eidors_objects.max_cache_size = limit;
      else
         retval= eidors_objects.max_cache_size;
      end

   case 'boost_priority'
      try
         retval= eidors_objects.cache_priority;
      catch
         retval= 0; % default priority
      end
      if nargin==2
         retval = retval + limit;
      end
      eidors_objects.cache_priority = retval;

   case 'show_objs'
      for i=1:length(times)
         fprintf('t=%9.7f b=%9.0d p=%02d, i=%03d: %s\n', rem(times(i),1), ... %today
             sizes(i), prios(i), priidx(i), objid{i} ); 
      end

   case 'clear_max'
% This will remove just in order of time.
%     remove_objids( objid, sizes,  find(cumsum(sizes) > limit) );
% Remove in order of time + priority
      tot= cumsum(sizes(priidx)); 
      tot= tot(priidx);
      remove_objids( objid, sizes,  find(tot > limit) );

   case 'clear_old'
      remove_objids( objid, sizes,  ...
        find(times < limit) );

   case 'clear_new'
      remove_objids( objid, sizes,  ...
        find(times > limit) );
   
   otherwise
      error('command %s not understood',command);
end

% priidx is priority index, where 1 prio is 0.1 of a day
function [objid, times, sizes, prios, priidx] = get_names_times;
   objid={}; times=[]; sizes=[]; pri=[]; prios=[]; priidx=[];
   global eidors_objects;
   if isempty(eidors_objects); return; end
   idx=1;
   for fn= fieldnames(eidors_objects)'
      fn1= fn{1};
      if fn1(1:3)== 'id_'
         objid{idx}= fn1;

         obj = getfield(eidors_objects, fn1 );
         times(idx) = obj.last_used;
         prios(idx) =   obj.priority;

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
   [jnk, priidx] = sort(-(times+0.1*prios));
   priidx(priidx)= 1:length(priidx);

function remove_objids( objid, sizes, idx)
   global eidors_objects;
   eidors_objects= rmfield( eidors_objects, objid( idx ) );
   eidors_msg('eidors_cache: removing %d objects with %d bytes', ...
          length(idx), sum(sizes(idx)), 3 );
