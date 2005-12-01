function Reg= ab_calc_tv_prior( inv_model );
% AB_CALC_TV_PRIOR calculate Total Variation image prior
% Reg= ab_calc_tv_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: ab_calc_tv_prior.m,v 1.2 2005-12-01 12:07:20 aadler Exp $

% Andrea's code requires a msh 
msh.TC = inv_model.fwd_model.elems';
msh.PC = inv_model.fwd_model.nodes';

dims= size(msh.PC,1);

if dims==2
    Reg= TV_operator_2D( msh );
elseif dims==3
    error('cant do this yet');
else
    error('problem dimension must be 2 or 3');
end
