function Reg= RtR_prior_elem2nodes( inv_model )
% RtR_PRIOR_ELEM2NODES: Convert elem to nodal RtR Image Prior
% Reg= RtR_prior_elem2nodes( inv_model)
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   inv_model.RtR_prior_elem2nodes.RtR_prior= func which calcs prior
%   inv_model.RtR_prior_elem2nodes.fwd_model= fwd model for elem2nodes

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

NtoE = mapper_nodes_elems( inv_model.RtR_prior_elem2nodes.fwd_model );

Reg_e= feval(inv_model.RtR_prior_elem2nodes.RtR_prior, inv_model);
Reg= NtoE*Reg_e*NtoE';

nn= size(NtoE,1);
% Sometimes NtoE is singular. We use a cheap test here, because it
%   is a structured matrix, and we don't want to use cond or condest
% Use a NOSER type prior, because its not quite clear how to 
%   regularize in this case.
singular = det(NtoE(:,1:nn))<1e-13;
if singular
   hp = 1e-6;
   dReg= hp*spdiags(Reg,0);
   Reg= Reg + spdiags(dReg,0,nn,nn);
end
