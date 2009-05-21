function J= filter_jacobian( fwd_model, img)
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

J0 = feval(fwd_model.filter_jacobian.jacobian, ...
           fwd_model, img);
F  = fwd_model.filter_jacobian.filter;

J= F*J0;

