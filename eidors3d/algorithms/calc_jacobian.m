function jacobian = calc_jacobian( fwd_model, img, varargin)
% CALC_JACOBIAN: calculate jacobian from a fwd_model object and an image
% 
%    jacobian= calc_jacobian( fwd_model, img, ...)
%
% jacobian  is a jacobian matrix
%      returned by the function fwd_model.solve
%
% fwd_model is a fwd_model structure
% img       is an image structure

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_jacobian.m,v 1.8 2007-08-29 09:15:30 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

jacobian = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.jacobian, img , varargin{:} );
