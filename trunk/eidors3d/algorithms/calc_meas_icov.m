function meas_icov = calc_meas_icov( inv_model )
% meas_icov = calc_meas_icov( inv_model )
% CALC_MEAS_ICOV: calculate inverse covariance of measurements
%   The meas_icov is a matrix n_meas x n_meas of the
%     inverse of measurement covariances. Normally measurements
%     are assumed to be independant, so the meas_icov is 
%     a diagonal matrix of 1/var(meas)
%
% calc_meas_icov can be called as
%    meas_icov= calc_meas_icov( inv_model )
%
% in each case it will call the inv_model.meas_icov
%  supplied by the user. If such a function is not available,
%  meas_icov will assume uniform independent measurement noise.
%
% meas_icov   is the calculated data prior
% inv_model    is an inv_model structure
%
% Many common EIT (and other regularized) algorithms do not
%  contain a meas_icov term. For these algorithms, this function
%  generates an approximation based on uniform noise.

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_meas_icov.m,v 1.17 2007-08-29 09:26:18 aadler Exp $

meas_icov = eidors_obj('get-cache', inv_model, 'meas_icov');

if ~isempty(meas_icov)
   eidors_msg('calc_meas_icov: using cached value', 3);
   return
end

if isfield(inv_model,'meas_icov')
   meas_icov= feval( inv_model.meas_icov, inv_model);
else
   meas_icov= default_meas_icov( inv_model );
end

eidors_obj('set-cache', inv_model, 'meas_icov', meas_icov);
eidors_msg('calc_meas_icov: setting cached value', 3);
 
% Calculate a data prior for an assumption of uniform noise
% on each channel
% 
function meas_icov = default_meas_icov( inv_model )

   fwd_model= inv_model.fwd_model;
   if isfield(fwd_model,'normalize_measurements')
      normalize = fwd_model.normalize_measurements;
   else
      normalize = 0;
   end

   n =  calc_n_meas( fwd_model );

   if ~normalize
      meas_icov= speye( n );
   else
      homg_data=  solve_homg_image( fwd_model );
      meas_icov = sparse(1:n, 1:n, ( 1./ homg_data.meas ).^2 );
   end

function n_meas = calc_n_meas( fwd_model )

   n_meas = 0;
   for i= 1:length(fwd_model.stimulation );
       n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern,1);
   end

% create homogeneous image + simulate data
function homg_data = solve_homg_image( fwd_mdl )
    n_elems= size( fwd_mdl.elems , 1);
    mat= ones( n_elems, 1);
    homg_img= eidors_obj('image', 'homogeneous image', ...
                         'elem_data', mat, 'fwd_model', fwd_mdl );
    homg_data=fwd_solve( homg_img);
