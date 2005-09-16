function jacobian = calc_jacobian( fwd_model, img)
% CALC_JACOBIAN: calculate jacobian from a fwd_model object and an image
% 
%    jacobian= calc_jacobian( fwd_model, img)
%
% it will call the fwd_model.solve
%
% jacobian  is a jacobian matrix
% fwd_model is a fwd_model structure
% img       is an image structure
%
% $Id: calc_jacobian.m,v 1.6 2005-09-16 03:17:46 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

jacobian = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.jacobian, img );
