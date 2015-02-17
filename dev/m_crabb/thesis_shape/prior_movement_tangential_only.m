function Reg= prior_movement_tangential_only( inv_model );
% Tikhonov movement prior with tangential only

pp= fwd_model_parameters( inv_model.fwd_model );

% calc movement portion
Reg= eye((pp.n_dims-1)*pp.n_elec);