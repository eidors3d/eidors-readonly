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
% $Id: calc_jacobian.m,v 1.4 2005-02-23 16:47:59 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

% TODO: caching jacobian depends on fwd_model and image
jacobian = eidors_obj('get-cache', fwd_model, 'jacobian');

if ~isempty(jacobian)
   eidors_msg('calc_jacobian: using cached value', 2);
   return
end

jacobian = feval( fwd_model.jacobian, fwd_model, img);
eidors_obj('set-cache', fwd_model, 'jacobian', jacobian);

eidors_msg('calc_jacobian: setting cached value', 2);
