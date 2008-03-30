function J = calc_jacobian( fwd_model, img)
% CALC_JACOBIAN: calculate jacobian from an inv_model
% 
%  J = calc_jacobian( fwd_model, img )
%  J = calc_jacobian( img )
%      calc Jacobian on fwd_model at conductivity given
%      in image (fwd_model is for forward and reconstruction)
%
% For reconstructions on dual meshes, the interpolation matrix
%    is defined as fwd_model.coarse2fine. This takes
%    coarse2fine * x_coarse = x_fine
%
% If the underlying jacobian calculator doesn't understand dual
%    meshes, then calc_jacobian will automatically postmultiply
%    by fwd_model.coarse2fine.
%
% img       is an image structure, with 'elem_data' or
%           'node_data' parameters

% (C) 2005-08 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_jacobian.m,v 1.22 2008-03-30 21:49:42 aadler Exp $

if nargin==1
   img       = fwd_model;
   fwd_model = img.fwd_model;
end

cache_obj= {fwd_model, img.elem_data};
J= eidors_obj('get-cache', cache_obj, 'jacobian');
if ~isempty(J)
   eidors_msg('calc_jacobian: using cached value', 3);
   return
end

J= feval(fwd_model.jacobian, fwd_model, img);

if isfield(fwd_model,'coarse2fine')
   c2f= fwd_model.coarse2fine;
   if size(J,2)==size(c2f,1)
%     calc_jacobian did not take into account the coarse2fine
      J=J*c2f;
   end
end

eidors_obj('set-cache', cache_obj, 'jacobian', J);
eidors_msg('calc_jacobian: setting cached value', 3);
