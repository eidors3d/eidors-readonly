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

warning('EIDORS:deprecated','BLD_MASTER_FULL is deprecated as of 07-Jun-2012. ');

dimen= size(vtx,2);
if dimen==2
   [Ef,D,Ela] = bld_master_full_2d(vtx,simp,mat,elec,zc);
elseif dimen==3
   [Ef,D,Ela] = bld_master_full_3d(vtx,simp,mat,elec,zc);
else
   error('not 2d or 3d');
end

function [Ef,D,Ela] = bld_master_full_2d(vtx,simp,mat,elec,zc);

[vr, vc] = size(vtx);
[sr, sc] = size(simp);
[er, ec] = size(elec);


if length(mat) ~= sr
   error('Invalid conductivity information for this mesh');
end


if length(zc) == er


%The column vector zc with the contact impedances in [Ohms] is required

[E,D,Ela] = bld_master(vtx,simp,mat);


E = full(E);

Ef = spalloc(vr+er,vr+er, er * vr);

Ef(1:vr,1:vr) = E;


while length(zc) ~= er
      disp(sprintf('Please enter the correct zc column vector with length: %d',er));
      %[zc] = contact_impedance;
end


for q=1:er
   
   tang_dist = 0;
   
   q_th_ele = elec(q,:);  % Select the row of nodes corresponding to the current electrode
   
   q_th_ele_zf = nonzeros(q_th_ele)'; % Extract the dummy "zero" nodal numbers
   
   for w=1:2:length(q_th_ele_zf)
      
      m = q_th_ele_zf(w);
      n = q_th_ele_zf(w+1);
      
      % This way m & n nodes belong to the edge tangent to the electrode and also at the same simplex.
      
      % We now evaluate the distance "tangential contact" between m & n 
      
      xm = vtx(m,1);
      ym = vtx(m,2); % m node coords
      xn = vtx(n,1);  
      yn = vtx(n,2); % n node coords
      
      [dist] = db2p(xm,ym,xn,yn); % distance mn
      
      cali_dist = dist ./ zc(q);  % coeficient for the distance mn
      
      tang_dist = tang_dist + cali_dist;
      
      % Start modifying "expanding" the E master matrix
      
      Ef(m,vr+q) = Ef(m,vr+q) - cali_dist/2 ; % Kv -> Ec  -> Vertical bar
      Ef(n,vr+q) = Ef(n,vr+q) - cali_dist/2 ; % Kv -> Ec
      
      Ef(vr+q,m) = Ef(vr+q,m) - cali_dist/2 ; % Kv' -> Ec' -> Horizontal bar
      Ef(vr+q,n) = Ef(vr+q,n) - cali_dist/2 ; % Kv' -> Ec'
      
      
      Ef(m,m) = Ef(m,m) + cali_dist/3; % Kz -> E -> Main bar
      Ef(n,n) = Ef(n,n) + cali_dist/3; % Kz -> E
      Ef(m,n) = Ef(m,n) + cali_dist/6; % Kz -> E
      Ef(n,m) = Ef(n,m) + cali_dist/6; % Kz -> E
      
   end % dealing with this electrode
   
   Ef(vr+q,vr+q) = Ef(vr+q,vr+q) + tang_dist;
   
end %for the whole set of electrodes

end

function [Ef,D,Ela] = bld_master_full_3d(vtx,simp,mat,elec,zc);
[vr, vc] = size(vtx);
[sr, sc] = size(simp);
[er, ec] = size(elec);


if length(mat) ~= sr
   error('Invalid conductivity information for this mesh');
end


[Ef,D,Ela] = bld_master(vtx,simp,mat);


% Add zeros so Ef is of size (vr+er) x (vr+er)
[Ef_i, Ef_j, Ef_s]= find( Ef );
Ef = sparse(Ef_i, Ef_j, Ef_s, vr+er, vr+er);


%Up to this point we have calculated the master matrix without the influence of contact impedance.

%The column vector zc with the contact 
%impedances in [Ohms] is required
if length(zc) ~= er
      error(sprintf('zc (=%d) should be equal to er (=%d)',length(zc),er));
end


for q=1:er
   
   tang_area = 0;
   
   q_th_ele = nonzeros(elec(q,:));  % Select the row of nodes corresponding to the current electrode
   
   if length(q_th_ele) ==1 % check if point electrode
      m = q_th_ele;
      cali_area = 1 / zc(q);
   
      tang_area = tang_area + cali_area;
      
      Ef(m,vr+q) = Ef(m,vr+q) - cali_area/2 ; 
      Ef(vr+q,m) = Ef(vr+q,m) - cali_area/2 ; 
      
      Ef(m,m) = Ef(m,m) + cali_area/2;

   else % not point electrode - use complete electrode model
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
   end % point electrode
   Ef(vr+q,vr+q) = Ef(vr+q,vr+q) + 0.5*tang_area;
   
end %for the whole set of electrodes

% calculate distance between two points
function [dist] = db2p(xa,ya,xb,yb);

   dist = sqrt((xb - xa).^2 + (yb - ya).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
