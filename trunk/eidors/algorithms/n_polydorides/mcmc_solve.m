function img= mcmc_solve( inv_model, data1, data2)
% NP_INV_SOLVE inverse solver for Nick Polydorides EIDORS3D code
% img= np_inv_solve( inv_model, data1, data2)
% img        => output image
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
% inv_model.parameters.max_iterations (default 1);
% inv_model.parameters.term_tolerance (default 1e-3);

% (C) 2007 Nick Polydorides. License: GPL version 2 or version 3
% $Id$

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

RtR = calc_RtR_prior( inv_model );
hp= calc_hyperparameter( inv_model );

img_bkgnd= calc_jacobian_bkgnd( inv_model );
J = calc_jacobian( inv_model.fwd_model, img_bkgnd);

sol= (J'*J +  hp^2*RtR)\(J' * dv );

%Set maximum number of iterations
try 
    max_iter= inv_model.parameters.max_iterations;
catch
    max_iter= 10;
end

likelihood= feval(inv_model.likelyhood_fcn, inv_model, sol, dv, J )



img.name= 'solved by mcmc_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;

