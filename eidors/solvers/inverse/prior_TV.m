function Reg= prior_TV( inv_model );
% PRIOR_TV calculate Total Variation image prior
% Reg= prior_TV( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

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

if isfield(inv_model.fwd_model,'coarse2fine')
     Reg = Reg*inv_model.fwd_model.coarse2fine;
end