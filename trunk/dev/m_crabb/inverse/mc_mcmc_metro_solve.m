function [img,cond_samp]= mc_mcmc_metro_solve( inv_model, sim_data, meas_data)
% MC_MCMC_METRO_SOLVE inverse solver 
% img= mc_mcmc_metro_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => simulated   data 
% data2      => measurement data
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

%Get parameters (default : 1 maxiter, 1e-3 tol,2 show_iter, 2 backtrack)
[n_samples,n_burn,n_show,one_step_reg,prop_var,noise_var,min_s,max_s,cond_frac] = get_parameters(inv_model);


%Calculate the background image and save as current image
img_bkgnd= calc_jacobian_bkgnd( inv_model ); img_cur=img_bkgnd.elem_data;
                     
%Calculate the prior matrix, hyperparameter 
R = ab_calc_tv_prior( inv_model );
hp= calc_hyperparameter( inv_model );

%Perform a one step GN reconstruction
J = calc_jacobian( inv_model.fwd_model, img_bkgnd);
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);
RtRtik=tikhonov_image_prior(inv_model);
grad_obj = J'*(-volt_diff_meas_sim);  
hess_obj = J'*J + one_step_reg^2*RtRtik;
grad_obj=-grad_obj; p_search = hess_obj \ grad_obj;
img_cur = img_cur + p_search; img_bkgnd.elem_data=img_cur;

%Show iterations and initialise an iteration
show_iter=1;
img_it.fwd_model= inv_model.fwd_model;
img_it.type='image';
%    img_it.elem_data=img_cur;
%    figure; show_fem(img_it,[1,0,0]);


%Initialise a solution matrix and put current image in first column
n_cond=size(img_cur,1); 
cond_samp=zeros(n_cond,n_samples+1); cond_samp(:,1)=img_cur;
n_cond_samp=ceil(n_cond*cond_frac); %no. pixel change per sample

%Prior
img_prior=img_cur;

    
%Caluclate the posterior
%Forward solve, calc prior and likelihood at first conductivity
img_bkgnd.elem_data=cond_samp(:,1);
sim_data=fwd_solve(img_bkgnd);
prior_b=calc_prior(R,cond_samp(:,1),hp);
likelihood_b=calc_likelihood(meas_data,sim_data,noise_var);        
%Posterior is product of prior and likelihood
post_b=prior_b*likelihood_b;

%Initialise the acceptance
acceptance=0;

%Loop through the samples
for i=1:n_samples           
    tic; fprintf('Sample %i of %i',i,n_samples);
    
    %Find a group of random conductivity pixels and then indices to change
    cond_rand=randperm(n_cond); 
    con_i=cond_rand(1:n_cond_samp); 
    
    %Cahce the old conductivity at the random pixels and perturb
    old_sigma=cond_samp(:,i);
    cond_samp(con_i,i)=cond_samp(con_i,i) + sqrt(prop_var)*randn(n_cond_samp,1);
               
    %Forward solve, calc prior and likelihood at new conductivity
    img_bkgnd.elem_data=cond_samp(:,i);
    sim_data=fwd_solve(img_bkgnd);
    prior_n=calc_prior(R,cond_samp(:,i),hp);    
    %Test for negative conducitivity
    for iii=1:n_cond_samp
        if( (cond_samp(con_i(iii),i) < min_s) || (cond_samp(con_i(iii),i) > max_s) )
        prior_n=0;
        end
    end
    likelihood_n=calc_likelihood(meas_data,sim_data,noise_var);
    %Posterior is product of prior and likelihood
    post_n=prior_n*likelihood_n;
        
    %3. See if new conductivity is accepted
    %Calculate the acceptance ratio and draw random uniform number
    alpha=min([1,post_n/post_b]); accept_rat=rand(1);
        
    if(alpha>=accept_rat) %Accept - Reassign posterior
        acceptance=acceptance+1;
        post_b=post_n;
    else %Reject - Resassign the conductivity
        cond_samp(:,i)=old_sigma;
    end      

    %Print the running acceptance rate
    fprintf(' has acceptance rate of %0.2f\n',acceptance/i); toc;

    %Update the sample
    cond_samp(:,i+1)=cond_samp(:,i);

    if(mod(i,n_show)==0 && show_iter==1)
        img_it.elem_data=mean(cond_samp(:,1:i),2);
        figure; show_fem(img_it,[1,0,0]);
    end
        
end

    
function prior= calc_prior(R_mat,img_val,alpha_param)
   %TV prior
%   TV_img_val_surf=R_mat*(img_val-img_prior);
%   TV_img_val=alpha_param*norm(TV_img_val_surf,1);
%   prior=exp(-TV_img_val);    
   
   %Tikhonov prior
   prior=exp(-alpha_param^2*norm(img_val-img_prior,2)^2);
   
end

function likelihood=calc_likelihood(meas_data,sim_data,noise_var)
    norm_diff=norm(meas_data.meas-sim_data.meas)^2;
    likelihood=exp(-0.5*norm_diff/noise_var);
    %Plot the prior to see proposal, set likelihood to 1
%    likelihood=1;
end
    

%Create a data structure to return
img.name= 'solved by mc_gibbs_solve';
img.fwd_model= inv_model.fwd_model;
img.type='image';

%Calculate the mean of samples (after burn)
img.elem_data=mean(cond_samp(:,n_burn:n_samples),2);




%Default parameters for the GN solver
function [n_sample, n_burn, n_show, one_step_reg, prop_var,noise_var,min_s,max_s,cond_frac] = get_parameters(inv_model)
   try
     n_sample= inv_model.parameters.n_sample;
   catch
     n_sample= 10^3;
   end

   try
     n_burn = inv_model.parameters.n_burn;
   catch
     n_burn = 200;
   end
   
   try
     n_show = inv_model.parameters.n_show;
   catch
     n_show = 50000;
   end
   
   try
     one_step_reg = inv_model.parameters.one_step_reg;
   catch
     one_step_reg = 10^-3;
   end     
   
   try
     prop_var = inv_model.parameters.prop_var;
   catch
     prop_var = 0.1;
   end     
      
   try
     noise_var = inv_model.parameters.noise_var;
   catch
     noise_var = 10^-3;
   end     
   
   try
     min_s = inv_model.parameters.min_s;
   catch
     min_s = 10^-3;
   end     
   
   try
     max_s = inv_model.parameters.max_s;
   catch
     max_s = 10^-3;
   end
   
   try
     cond_frac = inv_model.parameters.cond_frac;
   catch
     cond_frac = 0.01;
   end
   
   
end

end
