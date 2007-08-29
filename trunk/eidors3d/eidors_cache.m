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

% (C) 2005 Andy Adler. Licensed under GPL version 2
% $Id: eidors_cache.m,v 1.6 2007-08-29 09:10:13 aadler Exp $

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

[objid, times, sizes] = get_names_times;

switch command
   case {'clear_all','clear'}
      remove_objids( objid, sizes,  1:length(sizes) );

   case 'cache_size'
      if nargin==2
         eidors_objects.max_cache_size = limit;
      else
         retval= eidors_objects.max_cache_size;
      end

   case 'show_objs'
      for i=1:length(times)
         fprintf('t=%8.7f b=%9.0d: %s\n', rem(times(i),1), ... %today
             sizes(i), objid{i} ); 
      end

   case 'clear_max'
      remove_objids( objid, sizes,  ...
        find(cumsum(sizes) > limit) );

   case 'clear_old'
      remove_objids( objid, sizes,  ...
        find(cumsum(times) < limit) );

   case 'clear_new'
      remove_objids( objid, sizes,  ...
        find(cumsum(times) > limit) );
   
   otherwise
      error('command %s not understood',command);
end

function [objid, times, sizes] = get_names_times;
   objid={}; times=[]; sizes=[]; 
   global eidors_objects;
   if isempty(eidors_objects); return; end
   idx=1;
   for fn= fieldnames(eidors_objects)'
      fn1= fn{1};
      if fn1(1:3)== 'id_'
         objid{idx}= fn1;

         obj = getfield(eidors_objects, fn1 );
         times(idx) = obj.last_used;

         ww= whos('obj');
         sizes(idx) = ww.bytes;

         idx= idx+1;
      end
   end

   % sort from recent (high) to old (low)
   [times, idx] = sort(-times);
   objid= objid(idx);
   sizes= sizes(idx);

function remove_objids( objid, sizes, idx)
   global eidors_objects;
   eidors_objects= rmfield( eidors_objects, objid( idx ) );
   eidors_msg('eidors_cache: removing %d objects with %d bytes', ...
          length(idx), sum(sizes(idx)), 3 );
