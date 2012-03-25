function img = mc_GN_barrier_solve(inv_model,sim_data,meas_data)
%Do Gauss Netwon Method with barrier to ensure positvity of the con
%ductivity elements
%INPUT v_h,v_i - simulated/measured voltages
%      J - simulated Jacobian
%      tau - multipler of barrier parameter (strictly <1, default 0.1)
%      mu - barrier parameter
%      epsilon - term tolerance (about 1e-5-1e-3 - default 1e-3)
%      min_s - the minimal sigma per element (size(J,2))
%      max_s - the maximal   ""        ""        ""   
%      max_its - maximum iterations (default 1)
%
%% PARAMETERS - MAKE AUTOMATIC

%Convergence tolerance
[tol] = get_parameters(inv_model);

%Define the barrier parameters
tau=0.05; mu=1.0;

%Minimum and maximum allowed conductivity
min_s=0.8; max_s=2.0;

%Calculate the max iterations for Barrier method 
cnt=0; mu_it=mu; 
while(mu_it>1e-18); mu_it=mu_it*tau; cnt=cnt+1; end
maxiter=cnt;

%% DO THE PROGRAMMING

%Get number of conductivities and define upper and lower bound
n_cond=size(inv_model.fwd_model.elems,1);

%Calculate the background image and save as current image
img_bkgnd= calc_jacobian_bkgnd( inv_model );
img_bkgnd.elem_data=0.5*(max_s+min_s)*ones(n_cond,1);
img_cur=img_bkgnd.elem_data;
                     
%Calculate the Jacobian, prior matrix and hyperparameter
J = calc_jacobian( inv_model.fwd_model, img_bkgnd);
RtR = calc_RtR_prior( inv_model );
hp= calc_hyperparameter( inv_model );

%Calculate the voltage difference data (meas-sim)
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);

%Start the Gauss Newton iteration
for i=1:maxiter
     
       %Calculate the Barrier method matrices
       Psi=zeros(n_cond,1); A=zeros(n_cond,n_cond); Pi=zeros(n_cond,n_cond);
       for ii=1:n_cond          
           Psi(ii) = ( -2*img_cur(ii) + min_s + max_s )/ ...
              ( -img_cur(ii)^2 + ( min_s + max_s )*img_cur(ii) - min_s*max_s);
           
           Pi(ii,ii)=( mu/( -img_cur(ii)^2 + (min_s + max_s)*img_cur(ii) - min_s*max_s ))^2;
           
           A(ii,ii) = -2*img_cur(ii) + min_s + max_s;
       end
      
       %Gradient of objective function and negate for descent
       grad_beta_mu = J'*(-volt_diff_meas_sim) - mu*Psi;% + hp^2*RtR*(img_0-img_cur);
       grad_beta_mu=-grad_beta_mu;

       %Hessian of objective function
       hess_beta_mu = J'*J + hp^2*RtR +2*Pi + 1/mu*(A'*Pi*A);
       
       %Check condition number of Hessian
       cond_hess=cond(hess_beta_mu);
       if(cond_hess>1e6)
           %Calculate QR factorisation of A'
          [Q,R]=qr(A');
          
          %Calculate the modified gradient and negate
          grad_beta_mu=Q'*grad_beta_mu;
          
          %Calculate approximate Hessian for ill conditioning       
          G=1/mu*R*Pi*R'; hess_beta_mu=G;
          
          %Compute search direction and then multiply by orthogonal
          p_search_mu=hess_beta_mu\grad_beta_mu;
          p_search_mu=Q'\p_search_mu;    
       else                 
          %Compute the search direction - negate gradient and do search
          p_search_mu = hess_beta_mu \ grad_beta_mu; 
       end
        
       %% Backtracking line search
       %Line search parameters
       alpha=1.0; alpha_bls=0.1; beta_bls=0.5;
       
       %Create new candidate, forward solve and difference with measurements
       img_new = img_cur + alpha*p_search_mu;
       img_bkgnd.elem_data=img_new; 
       sim_data_new=fwd_solve(img_bkgnd.fwd_model,img_bkgnd);
       volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   

       %Calculate the functions for BLS
       beta_u_x_n_1 = beta_mu_f(volt_diff_meas_sim_new,img_new,min_s,max_s,mu);
       beta_u_x_n   = beta_mu_f(volt_diff_meas_sim,img_cur,min_s,max_s,mu);
       grad_u_x_n   = p_search_mu'*grad_beta_mu;
       while(beta_u_x_n_1 > beta_u_x_n + alpha_bls*alpha*grad_u_x_n)
           %Shrink the linesearch parameter and create new candidate
           alpha=alpha*beta_bls;
           
           %Create new candidate, forward solve and difference with measurements
           img_new = img_cur + alpha*p_search_mu;
       
           %Forward solve on new data and calc difference with measure
           img_bkgnd.elem_data=img_new; 
           sim_data_new=fwd_solve(img_bkgnd.fwd_model,img_bkgnd);
           volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
           
           %Calculate the functions for BLS
           beta_u_x_n_1 = beta_mu_f(volt_diff_meas_sim_new,img_new,min_s,max_s,mu);               
           if(alpha < 1e-16)
                error('Line Search failed')
           end
       end
           
       %Update the solution from the descent, decrease barrier and assign
       img_cur = img_new;
       img_bkgnd.elem_data=img_cur;
           
       %Resolve model, find difference data and test convergence
       sim_data_new=fwd_solve(img_bkgnd.fwd_model,img_bkgnd);
       volt_diff_meas_sim = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
       
       norm(volt_diff_meas_sim)      
       
       if norm(volt_diff_meas_sim)<tol; break; end  % test tolerance   
       
       %Recalculate the Jacobian and resolve forward problem
       J = calc_jacobian( inv_model.fwd_model, img_bkgnd);
       
       %Decrease the barrier parameter
       mu = mu*tau;
end

%Create a data structure to return
img.name= 'solved by mc_GN_barrier_solve';
img.elem_data = img_cur;
img.fwd_model= inv_model.fwd_model;
img.type='image';


%Objective function - voltage 2-norm (+ barrier??)
function [beta_mu]=beta_mu_f(diff_volt,img_data,min_s,max_s,mu)
    %Objective function
    beta_mu = 0.5*norm(diff_volt,2)^2;
    for iii=1:length(img_data)
       beta_mu=beta_mu - mu*log( -img_data(iii)^2+(min_s+max_s)*img_data(iii) - min_s*max_s);
    end

end


function [tol] = get_parameters(inv_model)
   try
     tol = inv_model.parameters.term_tolerance;
   catch
     tol= 1e-3;
   end 
end

end