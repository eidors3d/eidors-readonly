function jacobian = calc_jacobian( fwd_model, img)
% CALC_JACOBIAN: calculate jacobian from a fwd_model object and an image
% 
% calc_jacobian can be called as
%    jacobian= calc_jacobian( fwd_model, img)
% or
%    jacobian= calc_jacobian( img)
%
% in each case it will call the fwd_model.solve
%                        or img.fwd_model.solve method
%
% jacobian  is a jacobian matrix
% fwd_model is a fwd_model structure
% img       is an image structure
%
% $Id: calc_jacobian.m,v 1.2 2004-07-24 01:34:40 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

jacobian = eidors_obj('cache', img, 'jacobian');

if isempty(jacobian)
   jacobian = feval( fwd_model.jacobian, fwd_model, img);
   eidors_obj('cache', img, 'jacobian', jacobian);

   eidors_msg('calc_jacobian: setting cached value', 2);
else
   eidors_msg('calc_jacobian: using cached value', 2);
end
