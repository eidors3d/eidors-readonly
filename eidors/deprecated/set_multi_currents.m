function [I,Ib] = set_multi_currents(protocol,elec,vtx,gnd_ind,no_pl);
%function [I,Ib] = set_multi_currents(protocol,elec,vtx,gnd_ind,no_pl);
%
%This functions applies opposite or adjacent current patterns to each of
%the planes of the system simultaneously. 
%
%
%
%protocol= The selected protocol '{op}' or '{ad}'
%elec    = The electrodes
%vtx     = The vertices
%gnd_ind = the index of the ground node
%no_pl   = The number of planes
%Ib      = The current patterns
%I       = The RHS vectors, i.e., the current patterns padded with zeroes 
%          for the forward calculations

warning('EIDORS:deprecated','SET_MULTI_CURRENTS is deprecated as of 07-Jun-2012. ');

no_el = size(elec,1);

elpp = no_el/no_pl;
eld2 = elpp/2;


if protocol == '{op}'
   
d=eld2;
II = [];

	for j=1:no_pl
	Ib = [];
   
		for i=1:d
   	Ip = zeros(elpp,1);
	   Ip(i)= 1;
		Ip(i+eld2)= -1;
  		Ib = [Ib,Ip];
		end

II = [II;Ib];
end

I = zeros(size(vtx,1),size(Ib,2));
I = [I;II];

elseif protocol == '{ad}'
   
   d = elpp;
   
   II = [];
   
   for j=1:no_pl
       Ib = [];
      
       for i=1:d-1
          Ip = zeros(elpp,1);
          Ip(i)=1;
          Ip(i+1)=-1;
          Ib = [Ib,Ip];
       end
       
       lx = zeros(elpp,1);
       lx(end) = 1;
       lx(1) = -1;
       Ib = [Ib,lx];
       
   II = [II;Ib];
else
   error('protocol must be {ad} or {op}');
end

I = zeros(size(vtx,1),size(Ib,2));
I = [I;II];

end %protocol
       

Ib = I(size(vtx,1)+1:end,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
