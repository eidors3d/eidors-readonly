function [J] = jacobian_3d_fields(V,Ela,D,elec,vtx,simp,mat_ref,v_f,df);
% [J] = jacobian_3d_fields(V,elec,vtx,simp,mat_ref,v_f,df);
%
%calculates the Jacobian (sensitivity) matrix from V fwd
%
%J        = The Jacobian (sensitivity) matrix with respect to conductivity
%
%V        = forward solver voltage
%D,Ela    = parameters from the fem_master_full
%elec     = the electrodes matrix
%vtx      = The vertices matrix
%simp     = The simplices matrix
%mat_ref  = The reference conductivity vector
%v_f      = The measurement fields
%df       = Measurements per current pattern as used in v_f 

[vr,dim] = size(vtx);

el_no = size(elec,1);

if sum(df)~= size(v_f,2);
   error('Mismatched data input');
end

%Select the part referring to the interior nodes
V = V(1:vr,:);
v_f = v_f(1:vr,:);
   
J = zeros(sum(df),size(simp,1)); 
Jrow = zeros(1,size(simp,1));
cnt = 0;

el_idx= 1:dim:size(Ela,1); % Ela must be square
Vol_cond = diag(Ela(el_idx,el_idx)).* mat_ref;

   for p=1:size(V,2) 
     
      DV =  D*V(:,p); %Gradient of the current fields 
       
      for m=1:df(p) 
        
      
        Dvf = D*v_f(:,sum(df(1:p-1))+m); %Gradient of the measurement fields
              
        Jrow_x3 = Dvf .* DV ;  
        lJrow= length(Jrow_x3);
        Jrow_u = sum(reshape(Jrow_x3,dim,[]),1)'; % Works for 2D and 3D
        
%       Jrow_u = Jrow_x3(1:dim:lJrow) + Jrow_x3(2:dim:lJrow) + Jrow_x3(3:dim:lJrow);
%       Jrow = Jrow_u .* diag(Ela(1:3:size(Ela,1),1:3:size(Ela,2))); % Must account for background conductivity
        Jrow = Jrow_u .* Vol_cond;
        
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

