function Reg= prior_tikhonov( inv_model );
% PRIOR_TIKHONOV calculate image prior
% Reg= prior_tikhonov( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

%pp= fwd_model_parameters( inv_model.fwd_model );
% switch inv_model.type
%   case 'inv_model'; fwd_model = inv_model.fwd_model;
%   case 'fwd_model'; fwd_model = inv_model;
%   otherwise; error('PRIOR_TIKHONOV requires input type of inv_model or fwd_model');
% end

% if isfield( fwd_model, 'coarse2fine' )
%     no_dof = size(fwd_model.coarse2fine,2);
% else
%     no_dof = size(fwd_model.elems,1);
% end

img = physics_param_mapper(calc_jacobian_bkgnd(inv_model));
no_dof = numel(img.params);
Reg = speye( no_dof );

