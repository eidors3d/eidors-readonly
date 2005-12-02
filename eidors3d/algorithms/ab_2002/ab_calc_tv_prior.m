function Reg= ab_calc_tv_prior( inv_model );
% AB_CALC_TV_PRIOR calculate Total Variation image prior
% Reg= ab_calc_tv_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: ab_calc_tv_prior.m,v 1.4 2005-12-02 11:49:57 aadler Exp $

% Andrea's code requires a msh 
elem = inv_model.fwd_model.elems;
node = inv_model.fwd_model.nodes;

dims= size(node,2);

if dims==2
    msh.TC = elem';
    msh.PC = node';
    Reg= TV_operator_2D( msh );
elseif dims==3
    msh.elem_c = elem;
    msh.vtx_c  = node;
    Reg= TV_operator_3D( msh );
else
    error('problem dimension must be 2 or 3');
end
