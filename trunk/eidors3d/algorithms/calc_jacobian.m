function jacobian = calc_jacobian( fwd_model, img, varargin )
% CALC_JACOBIAN: calculate jacobian from an inv_model
% 
%  J = calc_jacobian( fwd_model, img )
%      calc Jacobian on fwd_model at conductivity given
%      in image (fwd_model is for forward and reconstruction)
%
%  J = calc_jacobian( fwd_model, img, rec_model )
%      calc Jacobian on fwd_model at conductivity given
%      in image (rec_model is for image reconstruction)
%
% fwd_model is a fwd_model structure
% rec_model is a fwd_model structure (but may not
%       have elems or electrodes)
% img       is an image structure, with 'elem_data' or
%           'node_data' parameters

% (C) 2005-08 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_jacobian.m,v 1.17 2008-03-13 20:42:11 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end

jacobian = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.jacobian, img , varargin{:} );
