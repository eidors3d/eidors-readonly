function Reg= tutorial210_cheat_tikhonov( inv_model )
% Reg= cheat_tikhonov( inv_model )
% Reg        => output regularization term
% Parameters:
%   elems    = inv_model.RtR_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.RtR_prior.cheat_weight;
%            => new weight to set elements to


pp= fwd_model_parameters( inv_model.fwd_model );
idx= 1:pp.n_elem;
weight= ones(1,pp.n_elem);
weight( inv_model.tutorial210_cheat_tikhonov.cheat_elements ) = ...
        inv_model.tutorial210_cheat_tikhonov.cheat_weight;

Reg = sparse( idx, idx, weight );
