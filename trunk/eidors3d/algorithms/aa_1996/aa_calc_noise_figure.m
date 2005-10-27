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
% $Id: aa_calc_noise_figure.m,v 1.13 2005-10-27 13:28:08 aadler Exp $

reqNF= inv_model.hyperparameter.noise_figure;

NFtable = eidors_obj('get-cache', inv_model, 'noise_figure_table');
if ~isempty(NFtable)
   % this would be sooo much easier if Matlab has assoc. arrays
   if any(NFtable(:,1) == reqNF)
       idx= find( NFtable(:,1) == reqNF);
       hparam= 10^NFtable( idx(1), 2);
       eidors_msg('aa_calc_noise_figure: using cached value', 2);
       return
   end
else
   NFtable= [];
end

startpoint = -5;
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
function [h_data, c_data, J]= simulate_targets( fwd_model, ctr_elems)

   %Step 1: homogeneous image
   sigma= ones( size(fwd_model.elems,1) ,1);

   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   h_data=fwd_solve( img );

   J = calc_jacobian( fwd_model, img);

   %Step 1: inhomogeneous image with contrast in centre
   delta = 1e-2;
   sigma(ctr_elems) = 1 + delta;
   img= eidors_obj('image', 'homogeneous image', ...
                   'elem_data', sigma, ...
                   'fwd_model', fwd_model );
   c_data=fwd_solve( img );

% calculate the noise figure for inv_model parameters
% based on the provided hyperparameter hp
function NF = calc_noise_figure( inv_model, hp)

   fwd_model= inv_model.fwd_model;
   pp= aa_fwd_parameters( fwd_model );

   [h_data, c_data, J]= simulate_targets( fwd_model, ...
        inv_model.hyperparameter.tgt_elems);
   if pp.normalize
      dva= 1 - c_data.meas ./ h_data.meas;
   else   
      dva= c_data.meas - h_data.meas;
   end

   R = calc_image_prior( inv_model );
   W = calc_data_prior( inv_model );

   % one step reconstruction matrix
   % for non-linear algorithms, we need the Taylor expansion at the soln
   RM= (J'*W*J +  hp*R)\J'*W;

%   NF = SNR_z / SNR_x
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(x0) / var(x) = sum(x0.^2) / trace(ABRnB'A')

   VOL2= pp.VOLUME'.^2;
   Rn= 1./diag(W); % assume independent noise, calc variance

   sig_data= sum( dva.^2);
   var_data= sum( Rn );
   sig_img = VOL2 * (RM * dva).^2;
   var_img = VOL2 * RM.^2 * Rn;

   NF = ( sig_data/ var_data ) / ( sig_img / var_img  );
   eidors_msg('calculating NF=%f hp=%g', NF, hp, 4);

   % For the record, the expression for var_img is derived as:
   % Equiv expresssions for var_img % given: A= diag(pp.VOLUME);
   % var_img= trace(A*RM*diag(Rn.^2)*RM'*A');
   % vv=A*RM*diag(Rn);var_img=trace(vv*vv'); var_img= sum(sum(vv.^2));
   % var_img= VOL2* (RM*diag(Rn)).^2
   % var_img= VOL2* RM.^2 * Rn.^2

   
