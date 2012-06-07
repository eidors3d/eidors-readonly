function [Er] = ref_master(E,vtx,gnd_ind,sch);
%function [Er] = ref_master(E,vtx,gnd_ind,sch);
%
%Applying reference to the system. Modifying the system matrix to
%preserve uniqueness.
%
%
%
%E       = The rank deficient by 1 system matrix
%Er      = The full rank matrix
%sch     = The grounding scheme:
%          0 for grounding node gnd_ind 
%          1 for grounding electrode gnd_ind
%gnd_ind = The ground index 

warning('EIDORS:deprecated','REF_MASTER is deprecated as of 07-Jun-2012. ');

[nv,jnk] = size(vtx);

if nargin < 4
   sch = 0;
end

if sch == 0 %Ground a surface node

Er = E;

Er(gnd_ind,:)= 0;					%zeros(1,mas_c);
Er(:,gnd_ind)= 0;					%zeros(mas_r,1);
Er(gnd_ind,gnd_ind) = 1;

else %Ground one of the boundary electrodes
   
Er = E;

if gnd_ind > size(E,1) - nv
    error('Grounding electrode index out of range')
end


Er(nv+gnd_ind,:)= 0;					%zeros(1,mas_c);
Er(:,nv+gnd_ind)= 0;					%zeros(mas_r,1);
Er(nv+gnd_ind,nv+gnd_ind) = 1;
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
