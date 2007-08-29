function NF = calc_noise_figure( inv_model, hp)
% AA_CALC_NOISE_FIGURE
% NF = calc_noise_figure( inv_model, hp)
% inv_model  => inverse model struct
%
% hp is specified, it will be used for the hyperparameter.
%    Otherwise the inv_model.hyperparameter will be used.
%
% Noise Figure must be defined for a specifc measurement
% In order to specify data, use
%     inv_model.hyperparameter.tgt_data.meas_t1
%     inv_model.hyperparameter.tgt_data.meas_t2
%   to use a temporal solver (or the Kalman filter), the
%   measurement to perform the NF calc must also be specified,
%   using:
%     inv_model.hyperparameter.tgt_data.meas_select
%   otherwise, the middle measurement will be used
%
% In order to automatically simulate data, specify tgt_elems,
%   containing a vector of elements to use
% 
%     inv_model.hyperparameter.tgt_elems
%

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_noise_figure.m,v 1.2 2007-08-29 09:10:10 aadler Exp $

% A 'proper' definition of noise power is:
%      NF = SNR_z / SNR_x
%    SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
%    SNR_x = sumsq(x0) / var(x) = sum(x0.^2) / trace(ABRnB'A')
% but his doesn't work, because the signal spreads like
% the amplitude, not like the power, thus
%      NF = SNR_z / SNR_x
%    SNR_z = sum(|z0|/len_z) / std(z/len_z)
%    SNR_x = sum(|x0|/len_x) / std(x/len_x)


   pp= aa_fwd_parameters( inv_model.fwd_model );

   if     isfield(inv_model.hyperparameter,'tgt_elems')
      [h_data, c_data]= simulate_targets( inv_model.fwd_model, ...
           inv_model.hyperparameter.tgt_elems);
   elseif isfield(inv_model.hyperparameter,'tgt_data')
      tgt_data= inv_model.hyperparameter.tgt_data;
      h_data= tgt_data.meas_t1;
      c_data= tgt_data.meas_t2;
   else
      error('unsure how to get data to measure signal');
   end

   % if hp is specified, then use that value
   if nargin>1
      inv_model= rmfield(inv_model,'hyperparameter');
      inv_model.hyperparameter.value= hp;
   else
      try
         hp = inv_model.hyperparameter.value;
      catch
         hp = NaN;
      end
   end
   if isfield(inv_model.hyperparameter,'func')
      if inv_model.hyperparameter.func == 'choose_noise_figure'
         error('specifying inv_model.hp.func == choose_noise_figure will recurse');
      end
   end

% To model std(z) we use z=z0+n
% so that std(z) = sqrt(var(z)) = sqrt(1/L * E[n'*n])
% we know a priori that the mean noise is zero, thus
% we do not need to divide by L-1 in the variance
% E[n'*n] = E[ trace(n*n') ] = trace( cov_N )
% Thus, std(z) = sqrt( trace( cov_N )/L ) 
%              = sqrt( mean( diag( cov_N )))
% And,  std(x) = sqrt( mean( diag( cov_X )))
%
% To model cov_N, we consider independant noise
%  on each channel, cov_N= N*N', N=diag( sigma_i )
% And,              cov_X= X*X', X=reconst(N)
% where X= reconst(mdl, z0,z0+N)
%
% To run efficiently mean(diag(cov_N))=mean(sum(N.^2,2))
% The sum over X is actually weighted by the volume of
%  each element, so VOL.^2*(sum(X.^2,2)

   % calculate signal
   d_len   = size(h_data,1);
   delta   = 1e-2* mean(h_data);
   c_noise = c_data*ones(1,d_len) + eye(d_len);
   h_full  = h_data*ones(1,d_len);

   sig_data = mean(abs( ...
         calc_difference_data( h_data, c_data , inv_model.fwd_model ) ...
                       ));
   var_data = mean(sum( ...
         calc_difference_data( h_full, c_noise, inv_model.fwd_model ) ...
                       .^2, 2)); 
                      

   % calculate image 
   % Note, this won't work if the algorithm output is not zero biased
   [img0, img0n] = get_images( inv_model, h_data, c_data, ...
                               h_full, c_noise);

   VOL = pp.VOLUME';
   sig_img= VOL*abs(img0.elem_data);
   var_img= VOL.^2*sum(img0n.elem_data.^2 ,2);
   
   NF = ( sig_data/ sqrt(var_data) ) / ( sig_img / sqrt(var_img)  );
   eidors_msg('calculating NF=%f hp=%g', NF, hp, 2);

   % For the record, the expression for var_img is derived as:
   % Equiv expresssions for var_img % given: A= diag(pp.VOLUME);
   % var_img= trace(A*RM*diag(Rn.^2)*RM'*A');
   % vv=A*RM*diag(Rn);var_img=trace(vv*vv'); var_img= sum(sum(vv.^2));
   % var_img= VOL2* (RM*diag(Rn)).^2
   % var_img= VOL2* RM.^2 * Rn.^2

   
% simulate homg data and a small target in centre
function [h_data, c_data]= simulate_targets( fwd_model, ctr_elems)

   homg= 1; % homogeneous conductivity level is 1

   %Step 1: homogeneous image
   sigma= homg*ones( size(fwd_model.elems,1) ,1);

   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   h_data=fwd_solve( img );
   h_data= h_data.meas;

   %Step 1: inhomogeneous image with contrast in centre
   delta = 1e-2;
   sigma(ctr_elems) = homg*(1 + delta);
   img.elem_data = sigma;
   c_data=fwd_solve( img );
   c_data= c_data.meas;

function [img0, img0n] = get_images( inv_model, h_data, c_data, ...
                               h_full, c_noise);
   if isa(inv_model.solve,'function_handle')
      solve= func2str(inv_model.solve);
   else
      solve= inv_model.solve;
   end

% Test for special functions and solve them specially
   switch solve
   case 'ab_tv_diff_solve'
      error('Dont know how to calculate TV noise figure')

   case 'inv_kalman_diff'
      inv_model.inv_kalman_diff.keep_K_k1= 1;
      stablize = 6;
      img0 = inv_solve( inv_model, h_data, ...
                                   c_data*ones(1,stablize) );
      K= img0.inv_kalman_diff.K_k1;
      img0.elem_data = K*calc_difference_data( h_data , c_data , inv_model.fwd_model);
      img0n.elem_data= K*calc_difference_data( h_full , c_noise, inv_model.fwd_model);

   otherwise
      img0 = inv_solve( inv_model, h_data, c_data);
      img0n= inv_solve( inv_model, h_full, c_noise);
   end
