function [Ef,D,Ela] = bld_master_full(vtx,simp,mat,elec,zc);
%function [Ef,D,Ela] = bld_master_full(vtx,simp,mat,elec,zc);
%
%System matrix assembling based on the complete electrode model. 
%This function is called within fem_master_full.
%
%
%
%Ef   = The UNreferenced system matrix.
%D    = The sgradients of the shape functions over each element.
%Ela  = Normalised volums of the elements
%vtx  = The vertices matrix. The coordinates of the nodes in 3D.
%simp = The simplices matrix. Unstructured tetrahedral.
%mat  = As for MATerial information. The conductivity vector.(isotropic)
%elec = The matrix that holds the boundary faces assigned as electrodes. Its typical
%       dimensions are [number of electrodes x 3*number of faces per electrode].
%zc   = The array of electrode contact impedances. 


[vr, vc] = size(vtx);
[sr, sc] = size(simp);
[er, ec] = size(elec);


if length(mat) ~= sr
   error('Invalid conductivity information for this mesh');
end


if length(zc) == er
%The column vector zc with the contact 
%impedances in [Ohms] is required


[Ef,D,Ela] = bld_master(vtx,simp,mat);


% Add zeros so Ef is of size (vr+er) x (vr+er)
[Ef_i, Ef_j, Ef_s]= find( Ef );
Ef = sparse(Ef_i, Ef_j, Ef_s, vr+er, vr+er);


%Up to this point we have calculated the master matrix without the influence of contact impedance.

if length(zc) ~= er
      error(sprintf('zc (=%d) should be equal to er (=%d)',length(zc),er));
end


for q=1:er
   
   tang_area = 0;
   
   q_th_ele = nonzeros(elec(q,:));  % Select the row of nodes corresponding to the current electrode
   
   for w=1:3:length(q_th_ele)
      
      m = q_th_ele(w);
      n = q_th_ele(w+1);
      l = q_th_ele(w+2);
      
        
      % This way m & n nodes belong to the edge tangential to the electrode 
      % and also at the same simplex.
      % We will now evaluate the distance "tangential contact area" between m,n & l 
      Are = triarea3d(vtx([m n l],:));
          
	% area mnl
      
      cali_area = (2*Are) ./ zc(q);  % coefficient for the area mnl
      %|J_k| = 2*Are  
      
      tang_area = tang_area + cali_area;
      
      % Start modifying "expanding" the E master matrix
      
      Ef(m,vr+q) = Ef(m,vr+q) - cali_area/6 ; % Kv -> Ec  -> Vertical bar
      Ef(n,vr+q) = Ef(n,vr+q) - cali_area/6 ; 
      Ef(l,vr+q) = Ef(l,vr+q) - cali_area/6 ;
            
      Ef(vr+q,m) = Ef(vr+q,m) - cali_area/6 ; % Kv' -> Ec' -> Horizontal bar
      Ef(vr+q,n) = Ef(vr+q,n) - cali_area/6 ; 
      Ef(vr+q,l) = Ef(vr+q,l) - cali_area/6 ;
      
      Ef(m,m) = Ef(m,m) + cali_area/12; % Kz -> E -> Main bar
      Ef(m,n) = Ef(m,n) + cali_area/24;       
      Ef(m,l) = Ef(m,l) + cali_area/24;
      
      Ef(n,m) = Ef(n,m) + cali_area/24;
      Ef(n,n) = Ef(n,n) + cali_area/12; 
      Ef(n,l) = Ef(n,l) + cali_area/24;
      
      Ef(l,m) = Ef(l,m) + cali_area/24;
      Ef(l,n) = Ef(l,n) + cali_area/24;
      Ef(l,l) = Ef(l,l) + cali_area/12;
    
      
   end % dealing with this electrode
   
   Ef(vr+q,vr+q) = Ef(vr+q,vr+q) + 0.5*tang_area;
   
end %for the whole set of electrodes

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
