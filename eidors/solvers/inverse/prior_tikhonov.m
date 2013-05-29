function Reg= prior_tikhonov( inv_model );
% PRIOR_TIKHONOV calculate image prior
% Reg= prior_tikhonov( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

%pp= fwd_model_parameters( inv_model.fwd_model );

if isfield( inv_model.fwd_model, 'coarse2fine' )
    no_dof = size(inv_model.fwd_model.coarse2fine,2);
else
    no_dof = size(inv_model.fwd_model.elems,1);
end

Reg = speye( no_dof );
