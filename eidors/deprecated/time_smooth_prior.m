function Reg= time_smooth_prior( inv_model );
% TIME_SMOOTH_PRIOR calculate image prior
% Reg= time_smooth_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% inv_model.time_smooth_prior.space_prior =
%          @space_prior_function
% inv_model.time_smooth_prior.time_steps  =
%          # of steps into future and past
% inv_model.time_smooth_prior.time_weight =  0..1
%    each step is weighted by time_weight^time_difference
%
% This image prior is intended to be used as
%  R'*R, but may be used as R for as well.
%
% The time smoothing prior penalizes non-smooth
% contributions in spatial and time directions
%
% The function of Reg is ||x-x_0||_Reg where 
% x is the image at 2*n+1 time slices concatenated
% vertially. x= [x_{j-n}; ... ; x_j ; ... x_{j+n} ]
%
% On a finite element mesh, we define the it as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','TIME_SMOOTH_PRIOR is deprecated as of 08-Jun-2012. Use PRIOR_TIME_SMOOTH instead.');
Reg = prior_time_smooth( inv_model);
