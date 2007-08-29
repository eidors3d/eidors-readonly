function img_bkgnd = calc_jacobian_bkgnd( inv_model )
% CALC_JACOBIAN_BKGND: calculate background image around
%    which initial estimate of jacobian is calculated
% 
% img_bkgnd = calc_jacobian_bkgnd( inv_model )
% inv_model   is an EIDORS fwd_model 
% img_bkgnd   is an EIDORS struct
%
% Usage: 
%      The background for calc_jacobian may be specified
%      as an estimated value
%  inv_model.jacobian_bkgnd.value;  % scalar OR
%                                   % vector of Nx1
%
% Usage:
%      The background may be calculated by a function
%  inv_model.jacobian_bkgnd.func;

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_jacobian_bkgnd.m,v 1.14 2007-08-29 09:20:55 aadler Exp $

img_bkgnd= eidors_obj('get-cache', inv_model, 'jacobian_bkgnd');
if ~isempty(img_bkgnd)
   eidors_msg('calc_jacobian_bkgnd: using cached value', 3);
   return
end

if isfield(inv_model.jacobian_bkgnd,'func')
   img_bkgnd= feval( inv_model.jacobian_bkgnd.func, inv_model );
else
   % allow bkgnd to be scalar or vector
   fwd_model= inv_model.fwd_model;
   bkgnd = ones(size(fwd_model.elems,1),1);
   bkgnd(:)= inv_model.jacobian_bkgnd.value;

   img_bkgnd= eidors_obj('image', 'homog image', ...
                         'elem_data', bkgnd, ...
                         'fwd_model', fwd_model );
end


eidors_obj('set-cache', inv_model, 'jacobian_bkgnd', img_bkgnd);
eidors_msg('jacobian_bkgnd: setting cached value', 3);
