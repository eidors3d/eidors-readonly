function Reg= time_smooth_prior( inv_model );
% TIME_SMOOTH_PRIOR calculate image prior
% Reg= time_smooth_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
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

% (C) 2006 Andy Adler. Licenced under the GPL Version 2
% $Id: time_smooth_prior.m,v 1.1 2006-07-27 11:33:41 aadler Exp $

pp= aa_fwd_parameters( inv_model.fwd_model );
ne = pp.n_elem;

% relative strengths of conductivity and movement priors
if ~isfield( inv_model,'time_smooth_prior')
   error('parameters must be specified for time_smooth_prior');
end

space_prior= inv_model.time_smooth_prior.space_prior;
time_weight= inv_model.time_smooth_prior.time_weight;
time_steps = inv_model.time_smooth_prior.time_steps;

space_Reg= feval(space_prior, inv_model);
time_Reg= -speye(ne);

Reg= [];
for i= -time_steps:time_steps;
   R_row= [];
   for j= -time_steps:time_steps;
      if i==j
         R_row= [R_row, space_Reg];
      else
         R_row= [R_row, time_Reg*time_weight*abs(i-j)];
      end
   end
   Reg= [Reg; R_row];
end
