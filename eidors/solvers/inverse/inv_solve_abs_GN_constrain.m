function [img,img_iteration] = inv_solve_abs_GN_constrain(inv_model,meas_data)
%Do Gauss Netwon Method with barrier to ensure positvity of the con
%ductivity elements
%INPUT v_h,v_i - simulated/measured voltages
%      inv_model - inverse model structure
%      tol - term tolerance (about 1e-5-1e-3 - default 1e-3)
%      min_s - the minimal sigma per element (size(J,2))
%      max_s - the maximal   ""        ""        ""   
%      max_its - maximum iterations (default 1)
%
%% PARAMETERS - MAKE AUTOMATIC
%
%M Crabb and N Polydorides- 29.06.2012
%TODO - Figure a nice interface to 
%       (i) Allow a best fitting homogeneous (user may not want this)
%       (ii) Allow a reference conductivity at each iteration
%       (iii) Allow a global reference conductivity
warning off backtrace
warning('EIDORS:experimental','%s is experimental, handle with care!',...
                upper('inv_solve_abs_GN_constrain'));
warning on backtrace

%Convergence tolerance
[maxiter, tol, min_s,max_s,rel_par,show_iter,bls,best_homog,best_homog_ref] = get_parameters(inv_model);

%Do not show iterations
if(show_iter==1); img_iteration=0; end

%Background image and simulate some data
img_bkgnd= calc_jacobian_bkgnd( inv_model );
img_cur=img_bkgnd.elem_data; 
sim_data=fwd_solve(img_bkgnd);

%"Best" fitting homogeneous?
if(best_homog==2) %BEST FIT
    sigma_opt=1/((sim_data.meas'*meas_data)/(sim_data.meas'*sim_data.meas));
    img_cur=img_cur*sigma_opt;
    img_bkgnd.elem_data=img_cur;
    sim_data=fwd_solve(img_bkgnd);
end
%TODO - Generalise this
img_ref=img_cur;

%Current image and the logistic equivalent
n_cond=length(img_cur);
log_img_cur=zeros(n_cond,1);
for i=1:n_cond
    log_img_cur(i)=inv_logistic_f(img_cur(i),min_s,max_s,rel_par);
end
                     
%Calculate the Jacobian, prior matrix and hyperparameter
RtR = calc_RtR_prior( inv_model );
hp= calc_hyperparameter( inv_model );

%Calculate the voltage difference data (meas-sim)
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);

%Start the Gauss Newton iteration
for i=1:maxiter
    %Print to screen if we want error
    if(show_iter==2);% && mod(i,ceil(maxiter/10))==0)       
        fprintf(1,'Error at iteration %i is %f\n',i,norm(volt_diff_meas_sim));
    end
    
    %Calculate the Jacobian
    J = calc_jacobian( img_bkgnd);
                
    %Compute the different vectors for method  (Polydorides 2012 Pg.10)
    d_s_d_m=zeros(n_cond,n_cond);
    for ii=1:length(log_img_cur)
        d_s_d_m(ii,ii)=d_logistic_f(log_img_cur(ii),min_s,max_s,rel_par);
    end
    
    %Multiply the partial derivative with Jacobian
    J=J*d_s_d_m;
    
    %Gradient of objective function (regularization term not needed)
    grad_obj = J'*(-volt_diff_meas_sim);
    
    %%TODO Implement with generic conductivity (see above)
    if(best_homog_ref==2)
        grad_obj=grad_obj - hp^2*RtR*(img_ref-img_cur);
    end

    %Hessian of objective function
    hess_obj = J'*J + hp^2*RtR;        
    
    %Compute search direction - negate gradient and do search    
    grad_obj=-grad_obj; p_search = hess_obj \ grad_obj;
    
    %% Backtracking line search??
    if(bls==2) %No linesearch
        %Update the constrained conductivity
        log_img_cur = log_img_cur + p_search; 
        %Change variables to normal conductivity
        for iii=1:n_cond
            img_cur(iii)=logistic_f(log_img_cur(iii),min_s,max_s,rel_par);
        end        
        img_bkgnd.elem_data=img_cur;
    else
       %Line search parameters
       alpha=1.0; alpha_bls=0.1; beta_bls=0.5;
       
       %Create new candidate, forward solve and difference with measurements
       log_img_new = log_img_cur + alpha*p_search; 
       
       %Change variables to normal conductivity
       for iii=1:n_cond
           img_new(iii)=logistic_f(log_img_new(iii),min_s,max_s,rel_par);
       end               
       img_bkgnd.elem_data=img_new'; 
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
           log_img_new = log_img_cur + alpha*p_search;
       
           %Change variables to normal conductivity
           for iii=1:n_cond
               img_new(iii)=logistic_f(log_img_new(iii),min_s,max_s,rel_par);
           end               
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
       log_img_cur=log_img_new; img_cur = img_new'; img_bkgnd.elem_data=img_cur;
    end
       
    %Resolve model, find difference data and test convergence
    sim_data_new=fwd_solve(img_bkgnd);
    volt_diff_meas_sim = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
    
    if norm(volt_diff_meas_sim)<tol; break; end  % test tolerance   
end

%Create a data structure to return
img.name= 'solved by mc_GN_box_solve';
img.elem_data = img_cur;
img.fwd_model= inv_model.fwd_model;
img.type='image';



%Objective function - voltage 2-norm (for linesearch)
function [beta]=beta_f(diff_volt)
    %Objective function
    beta = 0.5*norm(diff_volt,2)^2;
end

%Logistic function, its inverse and partial derivatives
function [logistic]=logistic_f(m_cur,min__s,max__s,relax_param)
    logistic = min__s + (max__s-min__s)/( 1 + exp(-m_cur/relax_param) );
end

function [d_logistic]=d_logistic_f(m_cur,min__s,max__s,relax_param)
    d_logistic = (max__s-min__s)/(   (1+exp(-m_cur/relax_param)) * ( (1+exp(m_cur/relax_param)) * relax_param));
end

function [inv_logistic]=inv_logistic_f(cond_cur,min__s,max__s,relax_param)
    inv_logistic = -relax_param*log( (cond_cur-max__s)/(min__s-cond_cur));
end


%Default parameters for the GN solver
function [maxiter, tol,min_s,max_s,rel_par,show_iter,bls,best_homog,best_homog_ref] = get_parameters(inv_model)
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
     min_s = inv_model.parameters.min_s;
   catch
     min_s= 1e-3;
   end
   
   try
     max_s = inv_model.parameters.max_s;
   catch
     tol= 1e3;
   end
   
   try
     rel_par = inv_model.parameters.rel_par;
   catch
     rel_par= 1.0;
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
