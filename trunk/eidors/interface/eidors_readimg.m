function img = eidors_readimg( fname, format )
% EIDORS readimg - read reconstructed image files from various EIT equipment
%    manufacturers
%
% Currently the list of supported file formats is:
%    - MCEIT (Goettingen / Viasys) "igt" file format 
%        format = "IGT" or "MCEIT"
%
% Usage
% [vv, auxdata ]= eit_readdata( fname, format )
%     vv      = measurements - data frames in each column
%     auxdata = auxillary data - if provided by system 
%     fname = file name
%
%  if format is unspecified, we attempt to autodetect

% (C) 2009 by Bartek Grychtol. Licensed under GPL v2 or v3
% $Id$ 


img.name = 'asdf'
img.type = 'image'
tempmdl = mk_common_gridmdl('b2d');
img.fwd_model = tempmdl.fwd_model;



loc.name = 'Dummy vector 1:856';
loc.elem_data = 1:856;
loc.fwd_model = img.fwd_model;
loc.type = 'image';