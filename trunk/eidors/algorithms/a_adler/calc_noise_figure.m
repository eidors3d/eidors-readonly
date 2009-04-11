function NF = calc_noise_figure( inv_model, hp, iterations)
% CALC_NOISE_FIGURE: calculate the noise amplification (NF) of an algorithm
% NF = calc_noise_figure( inv_model, hp, iterations)
%    inv_model  => inverse model object
%    hp         => value of hyperparameter to use (if not spec
%         then use the value of inv_model.hyperparameter.value)
%    iterations => number of iterations (default 10)
%       (for calculation of noise figure using random noise)
%
% hp is specified, it will be used for the hyperparameter.
%    Otherwise the inv_model.hyperparameter will be used.
%
% Noise Figure must be defined for a specific measurement
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
% $Id$

% A normal definition of noise power is based on power:
%      NF = SNR_z / SNR_x
%    SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
%    SNR_x = sumsq(x0) / var(x) = sum(x0.^2) / trace(ABRnB'A')
% but his doesn't work, because the signal spreads like
% the amplitude, not like the power, thus
%      NF = SNR_z / SNR_x
%    SNR_z = sum(|z0|/len_z) / std(z/len_z)
%    SNR_x = sum(|x0|/len_x) / std(x/len_x)


if nargin>=2
   inv_model.hyperparameter.value= hp;
   try; inv_model.hyperparameter = rmfield(inv_model.hyperparameter,'func'); end
else
   try hp= inv_model.hyperparameter.value; end
end
[inv_model, h_data, c_data, VOL] = process_parameters( inv_model );

%NF= nf_calc_use_matrix( inv_model, h_data, c_data, VOL);
%NF= nf_calc_iterate( inv_model, h_data, c_data, VOL); 
if nargin<3; N_RUNS= 10; end
NF= nf_calc_random( inv_model, h_data, c_data, VOL, N_RUNS);
eidors_msg('calculating NF=%f hp=%g', NF, hp, 2);

function [inv_model, h_data, c_data, VOL] = process_parameters( inv_model );
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

   if isfield(inv_model.hyperparameter,'func')
      funcname= inv_model.hyperparameter.func;
      if strcmp( class(funcname), 'function_handle')
         funcname= func2str(funcname);
      end

      if strcmp(funcname, 'choose_noise_figure')
         error('specifying inv_model.hp.func = choose_noise_figure will recurse');
      end
   end

   VOL = pp.VOLUME';

function NF= nf_calc_use_matrix( inv_model, h_data, c_data, VOL)
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

   i_len = length(img0.elem_data);
   sig_img= VOL*abs(img0.elem_data) / i_len;;
   var_img= VOL.^2*sum(img0n.elem_data.^2 ,2) / i_len;
   
   NF = ( sig_data/ sqrt(var_data) ) / ( sig_img / sqrt(var_img)  );

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
      if nargin>4
      img0n= inv_solve( inv_model, h_full, c_noise);
      end
   end

% OLD CODE - iterate
function NF= nf_calc_iterate( inv_model, h_data, c_data, VOL);
   % calculate signal
   d_len   = size(h_data,1);
   delta   = 1e-2* mean(h_data);
   sig_data = mean(abs( ...
         calc_difference_data( h_data, c_data , inv_model.fwd_model ) ...
                       ));
   % calculate image 
   % Note, this won't work if the algorithm output is not zero biased

   [img0] = get_images( inv_model, h_data, c_data);
%  sig_img= mean(VOL'.*abs(img0.elem_data));
   sig_img= VOL*abs(img0.elem_data) / length(img0.elem_data);

   % Now do noise
   var_data= 0;
   var_img = 0;
   for i=1:d_len
      this_noise = -ones(d_len, size(c_data,2))/(d_len-1);
      this_noise(i,:) = 1;
      c_noise = c_data + this_noise;
      [imgn] = get_images( inv_model, h_data, c_noise);
      if 1
         var_data = var_data + mean(sum( ...
            calc_difference_data( h_data, c_noise, inv_model.fwd_model ) ...
                          .^2, 2)); 
%        var_img= var_img +  mean( (VOL'.*imgn.elem_data).^2 ); 
         var_img= var_img +  (VOL.^2)*sum(imgn.elem_data.^2,2 ) / length(imgn.elem_data); 
      else
         % OLD APPROACH BASED ON variance, rather than matrix calcs
         var_data = var_data + var( ...
            calc_difference_data( h_data, c_noise, inv_model.fwd_model ) ...
                                 ); 
         var_img= var_img + var( VOL'.*imgn.elem_data ); 
      end
   end
   var_data = var_data / d_len;
   var_img  = var_img  / d_len;
   NF = ( sig_data/ sqrt(var_data) ) / ( sig_img / sqrt(var_img)  );

function NF= nf_calc_random( rec, vh, vi, VOL, N_RUNS);
   eidors_cache('boost_priority',-2); % low priority values

   imgr= inv_solve(rec, vh, vi);

   sig_ampl = mean( abs( VOL' .* imgr.elem_data )) / ...
              mean( abs( calc_difference_data( vh, vi, rec.fwd_model )));

% Estimate Signal Amplitude
   for i=1:N_RUNS
      vn= addnoise(vh, vi, 1.0);

      imgr= inv_solve(rec, vh, vn);
      noi_imag(i) = std( VOL' .* imgr.elem_data );
      noi_sgnl(i) = std( calc_difference_data( vh, vn, rec.fwd_model ));
   end
   noi_ampl = noi_imag./noi_sgnl;
   NF =  mean(noi_ampl/sig_ampl);
   SE =  std(noi_ampl/sig_ampl)/sqrt(N_RUNS);
   eidors_msg('NF= %f+/-%f\n', NF, SE, 1);

   eidors_cache('boost_priority',2);

   function noise= addnoise( vh, vi, SNR);
      if isstruct(vh); vh= vh.meas; end
      if isstruct(vi); vi= vi.meas; end
      noise = randn(size(vh));
      noise = noise*std(vh-vi)/std(noise);
      noise = vh + SNR*noise;

   
