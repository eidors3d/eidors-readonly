function [img,cond_samp]= mc_mcmc_abs_metro_solve( inv_model, meas_data)
% MC_MCMC_ABS_METRO_SOLVE inverse solver 
% img= mc_mcmc_abs_metro_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => simulated   data 
% data2      => measurement data
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

% STEP 1 - PARAMETERS + GLOBAL VARIABLES
%Parameters (overkill at moment)
[n_samples,n_burn,n_show,show_iter,one_step_reg,prop_var,noise_var,min_s,max_s,cond_frac] = get_parameters(inv_model);

%Background image, save as current image and simulate some data
img_bkgnd= calc_jacobian_bkgnd( inv_model ); img_cur=img_bkgnd.elem_data;
sim_data=fwd_solve(img_bkgnd);
                     
%Calculate prior matrix for conductivity
R = calc_R_prior( inv_model ); R_name=func2str(inv_model.R_prior);
hp= calc_hyperparameter( inv_model );

%STEP 2 - ONE STEP GN (faster convergence)
J = calc_jacobian( inv_model.fwd_model, img_bkgnd);
volt_diff_meas_sim = calc_difference_data( sim_data, meas_data, inv_model.fwd_model);
RtRtik=tikhonov_image_prior(inv_model);
grad_obj = J'*(-volt_diff_meas_sim);  
hess_obj = J'*J + one_step_reg^2*RtRtik;
grad_obj=-grad_obj; p_search = hess_obj \ grad_obj;
img_cur = img_cur + p_search; img_bkgnd.elem_data=img_cur;


%STEP 3 - Initialise solution matrix, current/prior/iteration images
n_cond=size(img_cur,1); 
if(memory_streamline==1) %Memory efficient method
    cond_samp=zeros(n_cond,2); 
else
    cond_samp=zeros(n_cond,n_samples+1); 
end
cond_samp(:,1)=img_cur;
n_cond_samp=ceil(n_cond*cond_frac); %no. pixel changes per sample
img_prior=img_cur; %Prior as first image? Should be input
acceptance=0; %acceptance rate
if(show_iter==1) 
    img_it.fwd_model= inv_model.fwd_model;
    img_it.type='image';
end

%STEP 4 - Calculate posterior at first step
%Forward solve, calc prior and likelihood at first conductivity
img_bkgnd.elem_data=cond_samp(:,1);
sim_data=fwd_solve(img_bkgnd);
prior_b=calc_prior(R,cond_samp(:,1),hp,R_name);
likelihood_b=calc_likelihood(meas_data,sim_data,noise_var);        
%Posterior is product of prior and likelihood
post_b=prior_b*likelihood_b;


%STEP 5 - Loop through samples, perturb and see if accepted
for i=1:n_samples           
    tic; fprintf('Sample %i of %i',i,n_samples);
    
    %Find a group of random conductivity pixels and then indices to change
    cond_rand=randperm(n_cond); 
    con_i=cond_rand(1:n_cond_samp); 
    
    %Cahce the old conductivity at the random pixels and perturb
    if(memory_streamline==1) %Memory efficient method
        old_sigma=cond_samp(:,1); cond_samp(con_i,1)=cond_samp(con_i,1) + sqrt(prop_var)*randn(n_cond_samp,1);
    else
        old_sigma=cond_samp(:,i); cond_samp(con_i,i)=cond_samp(con_i,i) + sqrt(prop_var)*randn(n_cond_samp,i);
    end       
    
    %Forward solve, calc prior and check if negative
    if(memory_streamline==1) %Memory efficient method
        img_bkgnd.elem_data=cond_samp(:,1); sim_data=fwd_solve(img_bkgnd); prior_n=calc_prior(R,cond_samp(:,1),hp,R_name);    
    else
        img_bkgnd.elem_data=cond_samp(:,i); sim_data=fwd_solve(img_bkgnd); prior_n=calc_prior(R,cond_samp(:,i),hp,R_name);    
    end     
     %Test for negative conducitivity
    for iii=1:n_cond_samp
        if(memory_streamline==1)
            if( (cond_samp(con_i(iii),1) < min_s) || (cond_samp(con_i(iii),1) > max_s) )
                prior_n=0;
            end
        else
            if( (cond_samp(con_i(iii),i) < min_s) || (cond_samp(con_i(iii),i) > max_s) )
                prior_n=0;
            end
        end
    end
    %Calc likelihood frome new data
    likelihood_n=calc_likelihood(meas_data,sim_data,noise_var);
    %Posterior is product of prior and likelihood
    post_n=prior_n*likelihood_n;
        
    %Determine if new conductivity accepted
    %Calculate acceptance ratio and draw random uniform number
    alpha=min([1,post_n/post_b]); accept_rat=rand(1);
        
    if(alpha>=accept_rat) %Accept - Reassign posterior
        acceptance=acceptance+1;
        post_b=post_n;
    else %Reject - Resassign the conductivity
        if(memory_streamline==1)
            cond_samp(:,1)=old_sigma;
        else
            cond_samp(:,i)=old_sigma;            
        end
    end      

    %Print acceptance rate every few samples
    fprintf(' has acceptance rate of %0.2f\n',acceptance/i); toc;
    
    %Update and show
    if(memory_streamline==1)    
        %After burn point start updating sample
        if(i>=n_burn)
            cond_samp(:,2)=cond_samp(:,2)+(cond_samp(:,1)-cond_samp(:,2))/(i+1-n_burn);
        end
        if(i >= n_burn && mod(i,n_show)==0 && show_iter==1)
            img_it.elem_data=cond_samp(:,2);
            figure; show_fem(img_it,[1,0,0]);
        end
    else
      %Update the sample and show 
        cond_samp(:,i+1)=cond_samp(:,i);
        if(mod(i,n_show)==0 && show_iter==1)
            img_it.elem_data=mean(cond_samp(:,1:i),2);
            figure; show_fem(img_it,[1,0,0]);
        end
    end 
end

    
function prior= calc_prior(R_mat,img_val,alpha_param,R_name) 
    %1 and 2 norm priors are fundamentally different
    if(strcmp(R_name,'prior_TV')) %1-norm type
        TV_img_val_surf=R_mat*(img_val-img_prior);
        TV_img_val=alpha_param*norm(TV_img_val_surf,1);
        prior=exp(-TV_img_val);    
    else %2-norm type
        prior=exp(-alpha_param^2*norm(R*(img_val-img_prior),2)^2);
    end
end

function likelihood=calc_likelihood(meas_data,sim_data,noise_var)
    %Assumed Gaussian distribution on voltages
    norm_diff=norm(meas_data.meas-sim_data.meas)^2;
    likelihood=exp(-0.5*norm_diff/noise_var);
end
    

%Create a data structure to return
img.name= 'solved by mc_gibbs_solve';
img.fwd_model= inv_model.fwd_model;
img.type='image';

%Calculate the mean of samples (after burn)
if(memory_streamline==1)
    img.elem_data=mean(cond_samp(:,2));
else
    img.elem_data=mean(cond_samp(:,n_burn:n_samples),2);
end


%Default parameters for the GN solver
function [n_sample, n_burn, n_show, show_iter, one_step_reg, prop_var,noise_var,min_s,max_s,cond_frac] = get_parameters(inv_model)
   try
     n_sample= inv_model.parameters.n_sample;
   catch
     n_sample= 100000;
   end

   try
     n_burn = inv_model.parameters.n_burn;
   catch
     n_burn = 10000;
   end
   
   try
     n_show = inv_model.parameters.n_show;
   catch
     n_show = 50000;
   end
   
   try
     show_iter = inv_model.parameters.show_iter;
   catch
     show_iter = 2;
   end
   
   try
     one_step_reg = inv_model.parameters.one_step_reg;
   catch
     one_step_reg = 10^-1;
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
     min_s = 10^-2;
   end     
   
   try
     max_s = inv_model.parameters.max_s;
   catch
     max_s = 10;
   end
   
   try
     cond_frac = inv_model.parameters.cond_frac;
   catch
     cond_frac = 0.01;
   end
   
   try 
     memory_streamline=inv_model.parameters.streamline;
   catch
     memory_streamline=1;    
   end

end

end
