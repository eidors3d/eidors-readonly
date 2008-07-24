function [J] = jacobian_3d_fields(V,Ela,D,elec,vtx,simp,mat_ref,v_f,df, c2f);
% [J] = jacobian_3d_fields(V,Ela,D,elec,vtx,simp,mat_ref,v_f,df, c2f);
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
%c2f      = Coarse to fine map between fine and coarse model (optional)

% $Id: jacobian_3d_fields.m,v 1.4 2008-05-09 22:38:00 aadler Exp $ 

[vr,dim] = size(vtx);

if sum(df)~= size(v_f,2);
   error('Mismatched data input');
end

%Select the part referring to the interior nodes
V = V(1:vr,:);
v_f = v_f(1:vr,:);

n_elem= size(simp,1);
%diag_Ela = diag(Ela(1:dim:size(Ela,1),1:dim:size(Ela,2)));
diag_Ela = diag(Ela(1:dim:size(Ela,1),1:dim:size(Ela,2)));
diag_Ela = spdiags(diag_Ela, 0, n_elem, n_elem);

if nargin>=10 % coarse2fine provided
   J = zeros(sum(df), size(c2f,2) );
   diag_Ela = diag_Ela*c2f;
else
   J = zeros(sum(df),n_elem );
end
cnt = 0;

for p=1:size(V,2)  % for each stimulation pattern

   DV =  D*V(:,p); %Gradient of the current fields 
   df_idx= sum(df(1:p-1));

   for m=1:df(p)   % for each measurement in this stim

     Dvf = D*v_f(:,df_idx + m); %Gradient of the measurement fields

     Jrow_x3 = Dvf .* DV ;  
     Jrow_u = sum(reshape(Jrow_x3,dim,[]),1); % Works for 2D and 3D now
     
%    Jrow = Jrow_u .* diag_Ela.';
     Jrow = Jrow_u * diag_Ela;
     
     cnt = cnt+1;
     J(cnt,:) = -Jrow;
     
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

