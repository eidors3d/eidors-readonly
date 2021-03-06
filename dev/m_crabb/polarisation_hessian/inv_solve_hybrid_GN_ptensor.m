function [ img,img_iteration ] = inv_solve_hybrid_GN_ptensor( inv_model, meas_data, hess_opts )
% INV_SOLVE_HYBRID_GN_PTENSOR Guass Newton inverse solver with added second
% derivatives on leading diagonal of Hessian, approximated by
% polarizability tensor expansion.
% img= mc_inv_solve_abs_GN( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => simulated   data 
% data2      => measurement data
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix
%M Crabb - 29.06.2012, F Watson 17.02.2017
%TODO - Figure a nice interface to 
%       (i) Allow a best fitting homogeneous (user may not want this)
%       (ii) Allow a reference conductivity at each iteration
%       (iii) Allow a global reference conductivity
warning off backtrace
warning('EIDORS:experimental','%s is experimental, handle with care!',...
                upper('inv_solve_abs_GN_prior'));
warning on backtrace


if nargin==2
    hess_opts = [];
end


if isfield(hess_opts, 'update_U0')
    update_U0 = hess_opts.update_U0;
else
    update_U0 = 0;
end
if isfield(hess_opts,'neumann')
    neumann = hesinput_argss_opts.neumann;
else
    neumann = 'freespace';
end


%Get parameters (default : 1 maxiter, 1e-3 tol,2 show_iter, 2 backtrack)
[maxiter, tol, show_iter,bls,best_homog,best_homog_ref] = get_parameters(inv_model);

%Do not show iterations
if(show_iter==1); img_iteration=0; end

%Background image, save as current image and simulate some data
img_bkgnd= calc_jacobian_bkgnd( inv_model );
img_cur=img_bkgnd.elem_data;
img_bkgnd.fwd_solve.get_all_meas = 1;
sim_data=fwd_solve(img_bkgnd);
if update_U0
    img_cur.fwd_solve.get_all_meas = 1;
end

%"Best" fitting homogeneous?
if(best_homog==2) %BEST FIT
    sigma_opt=1/((sim_data.meas'*meas_data)/(sim_data.meas'*sim_data.meas));
    img_cur=img_cur*sigma_opt;
    img_bkgnd.elem_data=img_cur;
    sim_data=fwd_solve(img_bkgnd);
end
%TODO - Generalise this
img_ref=img_cur;
                         
%Calculate the Jacobian, prior matrix and hyperparameter
RtR = calc_RtR_prior( inv_model );
hp= calc_hyperparameter( inv_model );

%Calculate the voltage difference data (meas-sim)
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);

% 
u0 = sim_data.volt;
DU0 = calc_grad_potential(img_bkgnd, u0);

%Start the Gauss Newton iteration
for i=1:maxiter
    %Print to screen if we want error
    if(show_iter==2);% && mod(i,ceil(maxiter/10)) == 0 )
        img_iteration{i}.error = norm(volt_diff_meas_sim);
        img_iteration{i}.name= sprintf('solved by mc_GN_solve iteration_%i',i);
        img_iteration{i}.elem_data = img_cur;
        img_iteration{i}.fwd_model= inv_model.fwd_model;
        img_iteration{i}.type='image';    
        
        fprintf(1,'Error at iteration %i is %f\n',i,norm(volt_diff_meas_sim));
%        figure;show_fem(img_iteration{i});
    end
    
    %Calculate the Jacobian
    J = calc_jacobian( img_bkgnd);

    %Gradient of objective function (regularization term not needed)
    %grad_obj = J'*W*(-volt_diff_meas_sim);
    grad_obj = J'*(-volt_diff_meas_sim);
    
    %%TODO Implement with generic conductivity (see above)
    if(best_homog_ref==2)
        grad_obj=grad_obj - hp^2*RtR*(img_ref-img_cur);
    end

    % ptensor
    if update_U0 && i>1
        u0 = sim_data.volt;
        DU0 = calc_grad_potential(img_bkgnd, u0);
    end
    [~,~,CC] = calc_phessian_obj(inv_model.fwd_model ,img_cur, DU0, volt_diff_meas_sim, neumann);
    
    %Hessian of objective function
    hess_obj = J'*J + hp^2*RtR + diag(CC);
    
    
    if ~exist('mod_chol','file')
        %Compute search direction - negate gradient and do search
        grad_obj=-grad_obj; p_search = hess_obj \ grad_obj;
        
    else
        % Modified Cholesky factorisation
        [L,D] = mod_chol(hess_obj);
        
        %Compute search direction - negate gradient and do search
        grad_obj=-grad_obj; p_search = L\D\(L.') \ grad_obj;
    end
    
    %% Backtracking line search??
    if(bls==2) %No linesearch
        img_cur = img_cur + p_search; img_bkgnd.elem_data=img_cur;
    else
       %Line search parameters
       alpha=1.0; alpha_bls=0.1; beta_bls=0.5;
       
       %Create new candidate, forward solve and difference with measurements
       img_new = img_cur + alpha*p_search; img_bkgnd.elem_data=img_new; 
       sim_data_new=fwd_solve(img_bkgnd);
       volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   

       %Calculate the functions for BLS
       beta_u_x_n_1 = beta_f(volt_diff_meas_sim_new);
       beta_u_x_n   = beta_f(volt_diff_meas_sim);
       grad_u_x_n   = p_search'*grad_obj;
       while(beta_u_x_n_1 > beta_u_x_n + alpha_bls*alpha*grad_u_x_n)
           %Shrink the linesearch parameter and create new candidate
           alpha=alpha*beta_bls;
           
           %Create new candidate, forward solve and difference with measurements
           img_new = img_cur + alpha*p_search;
       
           %Forward solve on new data and calc difference with measure
           img_bkgnd.elem_data=img_new; 
           sim_data_new=fwd_solve(img_bkgnd);
           volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
           
           %Calculate the functions for BLS
           beta_u_x_n_1 = beta_f(volt_diff_meas_sim_new);               
           if(alpha < 1e-16)
                error('Line Search failed');
           end
       end
       
       %Update the solution from the descent, decrease barrier and assign
       img_cur = img_new; img_bkgnd.elem_data=img_cur;
    end
       
    %Resolve model, find difference data and test convergence
    sim_data_new=fwd_solve(img_bkgnd);
    volt_diff_meas_sim = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
    
    if norm(volt_diff_meas_sim)<tol; break; end  % test tolerance   
end

%Create a data structure to return
img.name= 'solved by mc_GN_solve';
img.elem_data = img_cur;
img.fwd_model= inv_model.fwd_model;
img.type='image';



%Objective function - voltage 2-norm (for linesearch)
function [beta]=beta_f(diff_volt)
    %Objective function
    beta = 0.5*norm(diff_volt,2)^2;
end


%Default parameters for the GN solver
function [maxiter, tol,show_iter,bls,best_homog,best_homog_ref] = get_parameters(inv_model)
   try
     maxiter= inv_model.parameters.max_iterations;
   catch
     maxiter= 1;
   end

   try
     tol = inv_model.parameters.term_tolerance;
   catch
     tol= 1e-3;
   end
   
   try
     show_iter = inv_model.parameters.show_iterations;
   catch
     show_iter= 1;
   end
   
   try
      bls = inv_model.parameters.bls; 
   catch
      bls=1;
   end
   
   try
      best_homog = inv_model.parameters.best_homog; 
   catch
      best_homog=1;
   end
   
   try
      best_homog_ref = inv_model.parameters.best_homog_ref; 
   catch
      best_homog_ref=1;
   end
       
end

end
