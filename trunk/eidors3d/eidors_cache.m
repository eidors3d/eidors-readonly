function eidors_cache( command, varargin )
% Control eidors_caching
% Usage: eidors_cache( command, ... )
%
% USAGE:
%   eidors_cache clear all
%
%   eidors_cache clear all
%  
%
% (C) 2005 Andy Adler. Licensed under GPL version 2
% $Id: eidors_cache.m,v 1.3 2006-11-07 13:27:40 aadler Exp $

% Comments
% Want to clear specific structures
%      to clear old variables
%      to clear specific parts of structures

if nargin<1
    return;
end

global eidors_objects;
if strcmp( command , 'clear')
    eidors_objects= struct;
    clear eidors_objects;
else
    error('command not understood');
end
