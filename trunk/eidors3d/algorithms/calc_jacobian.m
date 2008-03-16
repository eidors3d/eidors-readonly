function J = calc_jacobian( fwd_model, img)
% CALC_JACOBIAN: calculate jacobian from an inv_model
% 
%  J = calc_jacobian( fwd_model, img )
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
% $Id: calc_jacobian.m,v 1.20 2008-03-16 11:07:44 aadler Exp $

if nargin>1
   img.fwd_model= fwd_model;
else
   img= fwd_model; % only img specified
end

J= eidors_obj('get-cache', img, 'jacobian');
if ~isempty(J)
   eidors_msg('calc_jacobian: using cached value', 3);
   return
end

J= feval(img.fwd_model.jacobian, img.fwd_model, img);

if isfield(img.fwd_model,'coarse2fine')
   c2f= img.fwd_model.coarse2fine;
   if size(J,2)==size(c2f,1)
%     calc_jacobian did not take into account the coarse2fine
      J=J*c2f;
   end
end

eidors_obj('set-cache', img, 'jacobian', J);
eidors_msg('calc_jacobian: setting cached value', 3);
