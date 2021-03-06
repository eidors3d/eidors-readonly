function J= filter_jacobian( varargin)
% FILTER_JACOBIAN: J= filter_jacobian( fwd_model, img)
%
% Filter a jacobian matrix by a specified filter function
% INPUT:
%  fwd_model = forward model
%  fwd_model.filter_jacobian.jacobian = actual jacobian function (J0)
%  fwd_model.filter_jacobian.filter   = Filter Matrix (F)
%  img = image background for jacobian calc
% OUTPUT:
%  J         = Jacobian matrix = F*J0

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','FILTER_JACOBIAN is deprecated as of 08-Jun-2012. Use JACOBIAN_FILTERED instead.');

if isfield(inv_model,'filter_jacobian');
  inv_model.jacobian_filtered = inv_model.filter_jacobian;
end

J = jacobian_filtered(varargin{:});

