function Reg= prior_movement_tangential( inv_model );
% PRIOR_MOVEMENT calculate image prior
% Reg= prior_movement( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   inv_model.image_prior.parameters(1) -> relative weighting
%     of movement vs image fraction of hyperparameter
%     => Default = 100
%   inv_model.prior_movement.RegC.func = Cond Reg fcn

% relative strengths of conductivity and movement priors
hp_move= inv_model.prior_movement.parameters(1);
pp= fwd_model_parameters( inv_model.fwd_model );

% calc conductivity portion
inv_model.RtR_prior = inv_model.prior_movement.RegC.func;
pp= fwd_model_parameters( inv_model.fwd_model );
RegC= calc_RtR_prior( inv_model); 
szC = size(RegC,1);
RegM = eye((pp.n_dims-1)*pp.n_elec);

%Create full regularisation
RegCM= sparse( szC, (pp.n_dims-1)*pp.n_elec );
Reg= [RegC,           RegCM;
      RegCM', hp_move^2*RegM ]; 