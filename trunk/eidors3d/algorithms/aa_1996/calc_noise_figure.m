function NF = calc_noise_figure( inv_model, hp)
% AA_CALC_NOISE_FIGURE
% NF = calc_noise_figure( inv_model, hp)
% inv_model  => inverse model struct
%
% hp is specified, it will be used for the hyperparameter.
%    Otherwise the inv_model.hyperparameter will be used.
%
% In order to use this function, it is necessary to specify
% inv_model.hyperparameter. has the following fields
% hpara.tgt_elems    = vector of element numbers of contrast in centre
%
% The NF parameter is modified from the definition in Adler & Guardo (1996).
%
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(A*x0) / var(A*x) = sum((A*x0).^2) / trace(ABRnB'A')
%   where Rn = noise covariance and A_ii = area of element i
% NF = SNR_z / SNR_x

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: calc_noise_figure.m,v 1.5 2006-11-24 04:21:27 aadler Exp $

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

   [h_data, c_data]= simulate_targets( inv_model.fwd_model, ...
        inv_model.hyperparameter.tgt_elems);

   % if hp is specified, then use that value
   if nargin>1
      inv_model= rmfield(inv_model,'hyperparameter');
      inv_model.hyperparameter.value= hp;
   else
      try
         hp = inv_model.hyperparameter.value;
      catch
         np = NaN;
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

   if inv_model.fwd_model.normalize_measurements
      sig_data= mean(abs(  c_data - h_data       ));
      var_data= mean(sum( (c_noise- h_full).^2, 2));
   else
      sig_data= mean(abs(  c_data ./ h_data -1       ));
      var_data= mean(sum( (c_noise./ h_full -1).^2, 2));
   end

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
      if inv_model.fwd_model.normalize_measurements
         img0.elem_data = K*( c_data  - h_data );
         img0n.elem_data= K*( c_noise - h_full );
      else
         img0.elem_data = K*( c_data ./ h_data - 1 );
         img0n.elem_data= K*( c_noise./ h_full - 1 );
      end

   otherwise
      img0 = inv_solve( inv_model, h_data, c_data);
      img0n= inv_solve( inv_model, h_full, c_noise);
   end
