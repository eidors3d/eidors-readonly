function J= perturb_jacobian( fwd_model, img)
% PERTURB_JACOBIAN: J= perturb_jacobian( fwd_model, img)
% Calculate Jacobian Matrix, based on small perturbations
%   in the forward model. This will tend to be slow, but
%   should be best used to 'sanity check' other code
%
% J         = Jacobian matrix
% fwd_model = forward model
% fwd_model.perturb_jacobian.delta   - delta perturbation to use
% img = image background for jacobian calc

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: perturb_jacobian.m,v 1.13 2007-08-30 03:37:02 aadler Exp $

if isfield(fwd_model,'perturb_jacobian')
   delta = fwd_model.perturb_jacobian.delta;
else
   delta= 1e-6; % tests indicate this is a good value
end

n_elem = size(fwd_model.elems,1);

% solve one time to get the size
d0= fwd_solve( img );
Jcol= perturb(img, 1, delta, d0);

J= zeros(length(Jcol), n_elem);
J(:,1)= Jcol;

for i=2:n_elem
  J(:,i)= perturb(img, i, delta, d0);
end

function Jcol= perturb( img, i, delta, d0)
   img.elem_data(i)= img.elem_data(i) + delta;
   di= fwd_solve( img );
   Jcol = (1/delta) * (di.meas - d0.meas);
