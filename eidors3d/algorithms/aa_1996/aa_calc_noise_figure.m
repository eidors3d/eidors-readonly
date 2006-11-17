function hparam= aa_calc_noise_figure( inv_model );
% AA_CALC_NOISE_FIGURE
% hparam= aa_calc_noise_figure( inv_model );
% inv_model  => inverse model struct
%
% In order to use this function, it is necessary to specify
% inv_model.hyperparameter. has the following fields
% hpara.func         = 'aa_calc_noise_figure';
% hpara.noise_figure = NF Value requested
% hpara.tgt_elems    = vector of element numbers of contrast in centre
%
% The NF parameter is modified from the definition in Adler & Guardo (1996).
% It is better to define it in terms of signal power, rather than amplitude.
%   measurements => z0 (Mx1), image elements => x0 (Nx1)
%
% given a (potentially nonlinear) reconstructor x=R(z), then
% the Taylor expansion x = R(z) ~= x0 + B*(z-z0) where B is NxM
%
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(A*x0) / var(A*x) = sum((A*x0).^2) / trace(ABRnB'A')
%   where Rn = noise covariance and A_ii = area of element i
% NF = SNR_z / SNR_x

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_calc_noise_figure.m,v 1.17 2006-11-17 14:31:50 aadler Exp $

reqNF= inv_model.hyperparameter.noise_figure;

NFtable = eidors_obj('get-cache', inv_model, 'noise_figure_table');
if ~isempty(NFtable)
   % this would be sooo much easier if Matlab has assoc. arrays
   idx= find( NFtable(:,1) == reqNF);
   if any(idx)
       hparam= 10^NFtable( idx(1), 2);
       eidors_msg('aa_calc_noise_figure: using cached value', 2);
       return
   end
else
   NFtable= [];
end

startpoint = -2;
opts = optimset('tolX',1e-4);
hparam= 10^fzero( @calc_log_NF, startpoint, opts, reqNF, inv_model );
   
NFtable = [NFtable; [reqNF, hparam] ];
eidors_obj('set-cache', inv_model, 'noise_figure_table', NFtable);
eidors_msg('aa_calc_noise_figure: setting cached value', 2);

% define a function that can be called by fzero. Also convert
% hparameter to log space to allow better searching by fzero
function out= calc_log_NF( log_hparam, reqNF, inv_model )
  out = calc_noise_figure( inv_model, 10^log_hparam ) - reqNF; 


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

% calculate the noise figure for inv_model parameters
% based on the provided hyperparameter hp
function NF = calc_noise_figure( inv_model, hp)

   pp= aa_fwd_parameters( inv_model.fwd_model );

   [h_data, c_data]= simulate_targets( inv_model.fwd_model, ...
        inv_model.hyperparameter.tgt_elems);


   inv_model= rmfield(inv_model,'hyperparameter');
   inv_model.hyperparameter.value= hp;

%   NF = SNR_z / SNR_x
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(x0) / var(x) = sum(x0.^2) / trace(ABRnB'A')

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

   
