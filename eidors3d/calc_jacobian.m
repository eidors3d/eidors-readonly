function jacobian = calc_jacobian( fwd_model, image)
% CALC_JACOBIAN: calculate jacobian from a fwd_model object and an image
% 
% calc_jacobian can be called as
%    jacobian= calc_jacobian( fwd_model, image)
% or
%    jacobian= calc_jacobian( image)
%
% in each case it will call the fwd_model.solve
%                      or image.fwd_model.solve method
%
% jacobian  is a jacobian matrix
% fwd_model is a fwd_model structure
% image     is an image structure
%
% $Id: calc_jacobian.m,v 1.1 2004-07-18 02:46:40 aadler Exp $

if nargin==1
   image= fwd_model;
   fwd_model= image.fwd_model;
end
jacobian= feval( fwd_model.jacobian, fwd_model, image);
