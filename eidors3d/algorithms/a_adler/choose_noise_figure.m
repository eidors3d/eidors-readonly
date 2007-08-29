function hparam= select_noise_figure( inv_model );
% SELECT_NOISE_FIGURE: select hyperparameter based on NF calculation
% hparam= select_noise_figure( inv_model );
% inv_model  => inverse model struct
%
% In order to use this function, it is necessary to specify
% inv_model.hyperparameter. has the following fields
% hpara.func         = @select_noise_figure;
% hpara.noise_figure = NF Value requested
% hpara.tgt_elems    = vector of element numbers of contrast in centre
%
% The NF parameter is from the definition in Adler & Guardo (1996).
%
% SNR_z = sumsq(z0) / var(z) = sum(z0.^2) / trace(Rn)
% SNR_x = sumsq(A*x0) / var(A*x) = sum((A*x0).^2) / trace(ABRnB'A')
%   where Rn = noise covariance and A_ii = area of element i
% NF = SNR_z / SNR_x

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: choose_noise_figure.m,v 1.2 2007-08-29 09:10:10 aadler Exp $

reqNF= inv_model.hyperparameter.noise_figure;

NFtable = eidors_obj('get-cache', inv_model, 'noise_figure_table');
if ~isempty(NFtable)
   % this would be sooo much easier if Matlab has assoc. arrays
   idx= find( NFtable(:,1) == reqNF);
   if any(idx)
       hparam= 10^NFtable( idx(1), 2);
       eidors_msg('select_noise_figure: using cached value', 2);
       return
   end
else
   NFtable= [];
end

startpoint = [4,-6]; % works better with a bracketed search
opts = optimset('tolX',1e-4);

% We don't want to cache any of these values
   pre_nf_timestamp = now;
hparam= 10^fzero( @calc_log_NF, startpoint, opts, reqNF, inv_model );
   eidors_cache('clear_new', pre_nf_timestamp);
   
NFtable = [NFtable; [reqNF, hparam] ];
eidors_obj('set-cache', inv_model, 'noise_figure_table', NFtable);
eidors_msg('select_noise_figure: setting cached value', 2);

% define a function that can be called by fzero. Also convert
% hparameter to log space to allow better searching by fzero
function out= calc_log_NF( log_hparam, reqNF, inv_model )
  out = calc_noise_figure( inv_model, 10^log_hparam ) - reqNF; 

