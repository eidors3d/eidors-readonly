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
% The NF parameter is defined in Adler & Guardo (1996), as
%   measurements => z (Mx1), image elements => x (Nx1)
%   NF = SNR_z / SNR_x
% SNR_z = mean(z) / std(z) = mean(z) /sqrt( M trace(Rn) )
% SNR_x = mean(x) / std(x) = mean(x) /sqrt( N trance(ABRnB'A)
%   where Rn = sigma_n x inv(W) = the noise covariance, iW= inv(W)
%     and A  = diag matrix s.t. A_ii = area of element i
%
% NF = mean(z)/mean(ABz) * sqrt(N/M) * sqrt( trace(ABiWB'A) / trance(iW) )
%
% given a reconstructor x=R(z), then
%   ABz = A R(z)
%
% NOTE: SNR _should_ be defined in terms of power! This defines
%       it in terms of images amplitude

% $Id: aa_calc_noise_figure.m,v 1.9 2005-09-13 02:37:54 aadler Exp $

% FIXME: this is a hack for now

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

   n_img = size(J,2);
   n_data= size(W,1);

   % one step reconstruction matrix
   RM= (J'*W*J +  hp*R)\J'*W;

   sig_data= sum(dva) / n_data;
   sig_img = ( pp.VOLUME' * RM * dva ) / n_img;

%  img= eidors_obj('image', 'x', ...
%                  'elem_data', RM*dva, 'fwd_model', fwd_model );
%  show_slices(img); pause
%  plot(RM*dva); pause

   var_data = trace(W)/n_data;
%  var_img  = sum( pp.VOLUME' * RM * DP ) / pp.n_elem;
% nf= sum(sig)*sqrt(sum(sum( ((AIRE*n_var').*Z).^2 ))) / ...
   var_img  = sum(sum( ((pp.VOLUME*diag(W)') .*RM).^2 )) / n_img;

% The NF parameter is calculated as follows
%   NF = SNR_z / SNR_x
% SNR_z = mean(z) / std(z) = mean(z) /sqrt( M trace(Rn) )
% SNR_x = mean(x) / std(x) = mean(x) /sqrt( N trance(ABRnB'A)

   NF = ( sig_data/sqrt(var_data) ) / ...
        ( sig_img /sqrt(var_img ) );
%  fprintf('%f - %f = %f\n', 1e6*sig_img, 1e6*var_img, NF );
   
