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

% (C) 2006 Andy Adler. Licenced under the GPL Version 2
% $Id: perturb_jacobian.m,v 1.1 2006-08-17 21:13:49 aadler Exp $

if isfield(fwd_model,'perturb_jacobian')
   delta = fwd_model.perturb_jacobian.delta;
else
   delta= 1e-6; % tests indicate this is a good value
end

pp= aa_fwd_parameters( fwd_model );

J = zeros( pp.n_meas, pp.n_elem );

elem_data = img.elem_data;
d0= fwd_solve( img );
for i=1:pp.n_elem
   img.elem_data   = elem_data;
   img.elem_data(i)= elem_data(i) + delta;
   di= fwd_solve( img );
   J(:,i) = (1/delta) * (d0.meas - di.meas);
end
