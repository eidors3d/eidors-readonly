function [img] = inv_solve_abs_GN_electrode_movement_sampling(inv_model,meas_data)
%Gauss Netwon Method with electrode positions and conductivity unknown with
%a line search to stop strange electrode perturbations 
%INPUT v_h,v_i - simulated/measured voltages
%      inv_model - inverse model structure
%      tol - term tolerance (about 1e-5-1e-3 - default 1e-3)
%      min_s - the minimal sigma per element (size(J,2))
%      max_s - the maximal   ""        ""        ""   
%      max_its - maximum iterations (default 1)
%

%Convergence tolerance
[maxiter, tol,alpha] = get_parameters(inv_model);

%Grab the forward model
fwd_model=inv_model.fwd_model;

%Calculate the useful electrode information i.e. tangents, normal etc.
elec_comp=calc_electrode_components(fwd_model);

%Calculate the DOFs for electrode movement (change basis where????)
dofelec=0;
for i=1:length(elec_comp)
   node_i_bound = length(elec_comp{i}.boundary_nodes);
   dofelec = dofelec+node_i_bound;
end

%Construct a Tikhonov matrix for this


%Get the number of element
n_elems=length(inv_model.fwd_model.elems(:,1));
n_elec=length(inv_model.fwd_model.electrode);
n_dim=length(inv_model.fwd_model.nodes(1,:));
for e=1:n_elec %Each electrode
    elec_nodes = inv_model.fwd_model.electrode(e).nodes;
    for i=1:n_dim %Each dimension
        pos(e,i) = mean(inv_model.fwd_model.nodes(elec_nodes,i),1);
    end
end

%Simulate some data on the background value (let this be reference atm)
img_bkgnd= calc_jacobian_bkgnd( inv_model );
img_cur=img_bkgnd.elem_data; 
sim_data=fwd_solve(img_bkgnd);

%Get the best fitting homogeneous conductivity
sigma_opt=1/((sim_data.meas'*meas_data)/(sim_data.meas'*sim_data.meas));
img_cur=img_cur*sigma_opt;
img_bkgnd.elem_data=img_cur;
sim_data=fwd_solve(img_bkgnd);

%Add on the electrode positions to the 'image'
for e=1:n_elec %Each electrode
    for i=1:n_dim %Each dimension
        img_cur(n_elems + (i-1)*n_elec + e)=pos(e,i);
    end
end
img_ref=img_cur;
                    
%Calculate the prior matrix - this has both the conductivity and movement
%as well as there relative weights embedded into it
RtR = calc_RtR_prior( inv_model );

%Calculate the overall regularisation parameter
hp= calc_hyperparameter( inv_model );

%Calculate the voltage difference data (meas-sim)
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);

%Start the Gauss Newton iteration
for i=1:maxiter
    %Print to screen if we want error
    fprintf(1,'Error at iteration %i is %f\n',i,beta_f(volt_diff_meas_sim,img_cur,img_ref,RtR));
    
    %Calculate the Jacobian with conductivity
    Jc = jacobian_adjoint( img_bkgnd); 
    
    %Calculate the Jacobian with movement
    Jm = jacobian_electrode_movement_sampling(img_bkgnd);
    
    %Construct the total Jacobian
    J = [Jc , Jm ];
    
    %Define a diagonal regularisation matrix
    reg_dim=size(J'*J,1); reg=eye(reg_dim); 
    
    
    %Gradient of objective function (regularization term not needed)
    grad_obj = J'*(-volt_diff_meas_sim) - hp^2*RtR*(img_ref-img_cur);

    %Hessian of objective function
    hess_obj = J'*J + hp^2*RtR;        
    
    %Compute search direction - negate gradient and do search    
    grad_obj=-grad_obj; p_search = hess_obj \ grad_obj;
    
    %% Backtracking line search
    %Line search parameters
    alpha=alpha; alpha_bls=0.1; beta_bls=0.5;
       
    %Create new candidate, forward solve and difference with measurements
    img_new = img_cur + alpha*p_search;    
    img_bkgnd.elem_data=img_new(1:n_elems);
    sim_data_new=fwd_solve(img_bkgnd);
    volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   

    %Armijo condition - fit data so just minimise the voltage part
    beta_u_x_n_1 = beta_f(volt_diff_meas_sim_new,img_new,img_ref,RtR); %Objective now just data
    beta_u_x_n   = beta_f(volt_diff_meas_sim,img_cur,img_ref,RtR);
    grad_u_x_n   = p_search'*grad_obj;
    while(beta_u_x_n_1 > beta_u_x_n + alpha_bls*alpha*grad_u_x_n)
        %Shrink the linesearch parameter and create new candidate
        alpha=alpha*beta_bls;
           
        %Create new candidate, forward solve and difference with measurements
        img_new = img_cur + alpha*p_search;       
        img_bkgnd.elem_data=img_new(1:n_elems); 
        sim_data_new=fwd_solve(img_bkgnd.fwd_model,img_bkgnd);
        volt_diff_meas_sim_new = calc_difference_data( sim_data_new, meas_data, inv_model.fwd_model);   
           
        %Calculate the functions for BLS
        beta_u_x_n_1 = beta_f(volt_diff_meas_sim_new,img_ref,img_new,RtR);               
        if(alpha < 1e-16)
            error('Line Search failed');
        end
    end
    
    fprintf(1,'BLS line search size at iteration %i is %f\n',i,alpha);
    
    %Update the solution from the descent, decrease barrier and assign
    img_cur = img_new; img_bkgnd.elem_data=img_cur(1:n_elems);
    
    img_show=img_cur; img_show.elem_data=img_bkgnd.elem_data;
    img_show.fwd_model= inv_model.fwd_model;
    img_show.type='image';    
    figure; show_fem(img_show,[1,0,0]);
       
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
function [beta]=beta_f(diff_volt,img_cur,img_ref,RtR)
    %Objective function
    beta = 0.5*norm(diff_volt,2)^2 + 0.5*(img_ref-img_cur)'*RtR*(img_ref-img_cur);
end

%Default parameters for the GN solver
function [maxiter,tol,alpha] = get_parameters(inv_model)
   try
     maxiter= inv_model.parameters.max_iterations;
   catch
     maxiter= 1;
   end
   
   try
     tol= inv_model.parameters.tolerance;
   catch
     tol= 10^-3;
   end
   
   try
     alpha=inv_model.parameters.alpha_bls;
   catch
     alpha=1;
   end
     
end

end
