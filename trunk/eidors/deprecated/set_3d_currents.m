function [I,Ib] = set_3d_currents(protocol,elec,vtx,gnd_ind,no_pl);
%function [I,Ib]=set_3d_currents(protocol,elec,vtx,gnd_ind,no_pl);
%
%This function sets current patterns in a system with (no_pl) planes of 
%equal number of electrodes according to "opposite" or "adjacent" protocols, 
%or their 3D similar.
%
%
%
%protocol= The selected protocol '{op}' or '{ad}'
%elec    = The electrodes (only the number of electrodes is used)
%vtx     = The vertices
%gnd_ind = the index of the ground node
%no_pl   = The number of planes
%Ib      = The current patterns
%I       = The RHS vectors, i.e., the current patterns padded with zeroes 


warning('EIDORS:deprecated','SET_3D_CURRENTS is deprecated as of 07-Jun-2012. ');

[vr,vc] = size(vtx);
   
[el_no,q] = size(elec); 

el_pp = el_no/no_pl;

a=1:el_no;

X = reshape(a,el_pp,no_pl)';

if protocol == '{op}'
   
   Ib = [];
   
	for i=1:no_pl
   
   	this_plane = X(i,:);
      
      for j=this_plane(1):this_plane(8)
         
         Ip = zeros(el_no,1);
         Ip(j) = 1;
         Ip(j+ el_pp/2) = -1;
         Ib = [Ib,Ip];
      end
   
   end 
   
    
   Is_supl = zeros(vr,size(Ib,2));

   I = [Is_supl;Ib];

   I(gnd_ind,:) = 0;

elseif protocol == '{ad}'
   
   Ib = [];
   
   for i=1:no_pl
   
     this_plane = X(i,:); 
   
        for j=this_plane(1):this_plane(el_pp-1)
           
           Ip = zeros(el_no,1);
           Ip(j) = 1;
           Ip(j+1) = -1;
           Ib =[Ib,Ip];
           
           if j==this_plane(el_pp-1) %the ring pattern
              
              Ip = zeros(el_no,1);
              
              Ip(j+1) = 1;
              Ip(this_plane(1)) = -1;
              Ib = [Ib,Ip];
           end
           
           
        end
        
     end
     
    Is_supl = zeros(vr,size(Ib,2));

    I = [Is_supl;Ib];

    I(gnd_ind,:) = 0;

else
   error(['protocol ',protocol,' is not recognized']);
end %protocol
         


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
