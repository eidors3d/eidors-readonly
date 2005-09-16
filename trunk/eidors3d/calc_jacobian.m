function jacobian = calc_jacobian( fwd_model, img, cachename)
% CALC_JACOBIAN: calculate jacobian from a fwd_model object and an image
% 
%    jacobian= calc_jacobian( fwd_model, img, cachename)
%
%
% jacobian  is a jacobian matrix
%      returned by the function fwd_model.solve
%
% fwd_model is a fwd_model structure
% img       is an image structure
% cachename is an explicit file to cache to 'cache:text' [OPTIONAL]
%
% $Id: calc_jacobian.m,v 1.7 2005-09-16 03:36:50 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

jacobian = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.jacobian, img );
