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
% $Id: calc_noise_figure.m,v 1.1 2006-11-17 14:47:38 aadler Exp $

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
   end
   if exist(inv_model.hyperparameter,'func')
      if inv_model.hyperparameter.func == 'choose_noise_figure'
         error('specifying inv_model.hp.func == choose_noise_figure will recurse');
      end
   end

   % calculate signal
   d_len   = size(h_data,1);
   delta   = 1e-2* mean(h_data);
   c_noise = c_data*ones(1,d_len) + eye(d_len);
   h_full  = h_data*ones(1,d_len);
   if inv_model.fwd_model.normalize_measurements
      sig_data= norm(c_data - h_data);
      var_data = trace((c_noise - h_full).^2);
   else
      sig_data= norm( c_data ./ h_data - 1);
      var_data = trace((c_noise./ h_full - 1).^2);
   end

   % calculate image 
   % Note, this won't work if the algorithm output is not zero biased
      img0 = inv_solve( inv_model, h_data, c_data);
      img0n= inv_solve( inv_model, h_full, c_noise);
   if 0
      VOL2= pp.VOLUME'.^2;
      sig_img = VOL2*img0.elem_data.^2;
      var_img = sum(VOL2*img0n.elem_data.^2); % trace weighted by area
   else
      if inv_model.fwd_model.normalize_measurements
         sig_data= mean(abs(c_data - h_data));
      else
         sig_data= mean(abs( c_data ./ h_data - 1));
      end
      var_data= sqrt(var_data);

      VOL = pp.VOLUME';
      sig_img= VOL*img0.elem_data;
      var_img = sqrt(sum(VOL.^2*img0n.elem_data.^2));
   end
   
   NF = ( sig_data/ var_data ) / ( sig_img / var_img  );
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
