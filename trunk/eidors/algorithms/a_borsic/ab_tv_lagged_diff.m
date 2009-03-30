function img=ab_tv_lagged_diff( inv_model, data1, data2)
% AB_TV_LAGGED_DIFF Lagged Diffusivity Inverse Solver
% img= ab_inv_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% Parameters:
%  max_iters =  inv_model.parameters.max_iteration (default 10)
%      Max number of iterations before stopping
%  min change = inv_model.parameters.min_change   (default 0)
%      Min Change in objective fcn (norm(y-Jx)^2 + hp*TV(x)) before stopping
%  beta      =  inv_model.ab_tv_lagged_diff.beta   (default 1e-3)
% beta is the parameter that smooths the TV functional

% (C) 2008 Andrea Borsic. License: GPL version 2 or version 3
% $Id%

try    max_iter = inv_model.parameters.max_iterations;
catch  max_iter = 10;
end

try    min_change = inv_model.parameters.min_change;
catch  min_change = 0;
end

try    beta = inv_model.ab_tv_lagged_diff.beta; 
catch  beta = 1e-3;
end


fwd_model= inv_model.fwd_model;
d=calc_difference_data( data1, data2, inv_model.fwd_model);

L=calc_R_prior( inv_model );

img_bkgnd=calc_jacobian_bkgnd( inv_model );
J=calc_jacobian( fwd_model, img_bkgnd);

alpha=calc_hyperparameter( inv_model );


delta_sigma = zeros(size(J,2),1); % we start from no initial difference

Obj_Fcn = inf; %initial value

for k=1:max_iter
 
    dv =  J*delta_sigma - d;

% STOP IF OBJECTIVE FCN IS NOT CHANGING by more than min_change
    Obj_Fcnk = norm(dv)^2 + alpha*sum(abs( delta_sigma ));
    delta_ObjFcn = abs(Obj_Fcnk/Obj_Fcn - 1);
%   fprintf('%d %g %g %g\n',k, Obj_Fcnk, Obj_Fcn, delta_ObjFcn);
    if delta_ObjFcn < min_change; 
       eidors_msg('Lagged_diff: Breaking at iteration %d',k,2);
       break;
    end
    Obj_Fcn = Obj_Fcnk;
    
    E= sqrt((L*delta_sigma).^2+beta);
    inv_E= spdiags( 1./E, 0, length(E), length(E));

    phi1=J'*dv+ alpha*L'*inv_E*L*delta_sigma;
    phi2=J'*J + alpha*L'*inv_E*L;

    upd=-(phi2)\phi1;
    
    delta_sigma=delta_sigma+upd;

    
end % for k

% create a data structure to return
img.name = 'ab_tv_lagged_diff';
img.elem_data = delta_sigma;
img.fwd_model = fwd_model;
