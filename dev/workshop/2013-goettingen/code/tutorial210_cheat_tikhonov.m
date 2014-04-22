function Reg= tutorial210_cheat_tikhonov( inv_model )
% Reg= cheat_tikhonov( inv_model )
% Reg        => output regularization term
% Parameters:
%   elems    = inv_model.RtR_prior.cheat_elements;
%            => elements weights to modify
%   weight   = inv_model.RtR_prior.cheat_weight;
%            => new weight to set elements to


ne = num_elems(inv_model.fwd_model);

weight= ones(1,ne);
weight( inv_model.tutorial210_cheat_tikhonov.cheat_elements ) = ...
        inv_model.tutorial210_cheat_tikhonov.cheat_weight;

Reg = spdiag( weight );
