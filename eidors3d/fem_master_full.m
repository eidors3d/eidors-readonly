function [E,D,Ela,pp] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,sym);
%function [E,D,Ela,pp] = fem_master_full(vtx,simp,mat,gnd_ind,elec,zc,sym);
%
%Builds up the system matrix based on the complete electrode model. E is not 
%yet permuted. To permute E -> E(pp,pp) as in forward_solver.
%
%
%
%E       = The full rank system matrix based on the 3D complete electrode model
%D       = The sgradients of the shape functions over each element.
%Ela     = Normalised volums of the elements
%pp      = Column permutation vector, for more help type help symmmd 
%vtx     = The vertices matrix
%simp    = The simplices matrix
%mat     = The conductivity vector
%gnd_ind = The index of the ground node
%elec    = The bounary electrodes matrix
%zc      = The contact impedance vector, satisfying size(elec,1) = length(zc)
%sym     = Column permutation of E, either '{y}' to opt or '{n}' to avoid.       

     
   [Ef,D,Ela] = bld_master_full(vtx,simp,mat,elec,zc); 
   
   [E] = ref_master(Ef,vtx,gnd_ind);  
   
 if sym == '{n}'
    pp = 1:size(E,1);
end
 if sym == '{y}'
    pp = symmmd(E);
end
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

