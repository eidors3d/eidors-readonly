function [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,df,tol,sym);
%function [J] = jacobian_3d(I,elec,vtx,simp,gnd_ind,mat_ref,zc,v_f,df,tol,sym);
%
%This function calculates the Jacobian (sensitivity) matrix of the system at
%mat_ref. 
%
%
%
%I        = The currents used
%elec     = the electrodes matrix
%vtx      = The vertices matrix
%simp     = The simplices matrix
%gnd_ind  = The ground index (node)
%mat_ref  = The reference conductivity vector
%zc       = The electrode contact impedance vector
%IntGrad  = The integrals of the gradients
%v_f      = The measurement fields
%df       = Measurements per current pattern as used in v_f 
%tol      = Tolerance 
%J        = The Jacobian (sensitivity) matrix with respect to conductivity

[vr,vc] = size(vtx);
[sr,sc] = size(simp);

el_no = size(elec,1);

if sum(df)~= size(v_f,2);
   error('Mismatched data input');
end

[E,D,Ela,pp] = fem_master_full(vtx,simp,mat_ref,gnd_ind,elec,zc,sym);

[V] = forward_solver(vtx,E,I,tol,pp);

%Select the part referring to the interior nodes
V = V(1:vr,:);
v_f = v_f(1:vr,:);
   
J = zeros(sum(df),size(simp,1)); 
Jrow = zeros(1,size(simp,1));
cnt = 0;

   for p=1:size(V,2) 
     
      DV =  D*V(:,p); %Gradient of the current fields 
       
      for m=1:df(p) 
        
      
        Dvf = D*v_f(:,sum(df(1:p-1))+m); %Gradient of the measurement fields
              
        Jrow_x3 = Dvf .* DV ;  
        Jrow_u = Jrow_x3(1:3:end) + Jrow_x3(2:3:end) + Jrow_x3(3:3:end);
        
        Jrow = Jrow_u .* diag(Ela(1:3:end,1:3:end));
        
        cnt = cnt+1;
        J(cnt,:) = -Jrow.';
        Jrow = zeros(1,size(simp,1));
        
    end %m
       
   end %p
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides and W.R.B. Lionheart 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%