function y = convert_units(varargin)
%CONVERT_IMG_UNITS change image data units 
%  img = convert_img_units(img,new_unit) converts img.elem_data or img.node_data
%  expressed in img.current_params into different units for the same 
%  physical property
% 
% Examples: 
%  img = mk_image(mdl, 2, 'resisitivity');
%  img = convert_img_units(img, 'conductivity');
%
%  img = mk_image(mdl, 1);
%  img.current_params = 'log_resistivity';
%  img = convert_img_units(img, 'conductivity');

% (C) 2012-2014 Alistair Boyle and Bartlomiej Grychtol
% License: GPL version 2 or 3
% $Id$

warning('EIDORS:deprecated','CONVERT_UNITS is deprecated as of 30-Apr-2014. Use CONVERT_IMG_UNITS instead.');
y = convert_img_units(varargin{:});
