function [IntGrad] = integrofgrad(vtx,simp,mat_ref);
%function [IntGrad] = integrofgrad(vtx,simp,mat_ref);
%
%function that calculates the integral of the gradients for first order
%tetrahedral elements. Required for the calculation of the Jacobian.
%
%
%
%vtx     = The vertices matrix
%simp    = The simplices matrix
%mat_ref = The reference conductivity vector
%IntGrad = The intgral of the gradients



warning('EIDORS:deprecated','INTEGROFGRAD is deprecated as of 06-Jun-2012. ');

[vr,vc] = size(vtx);
[sr,sc] = size(simp);

if sr ~= length(mat_ref)
   error('Mismatched data entered')
end

if sc ~= 4
   error('Only first order tetrahedral elements supported');
end


%w = [1/24*ones(4,1)];

IntGrad = sparse(vr^2,sr);

for i=1:size(simp,1)

   vv = [];

   for j=1:size(simp,2)
      vv = [vv;vtx(simp(i,j),:)];
   end

   interp_fun = [-1 1 0 0; -1 0 1 0; -1 0 0 1];

   Jts = interp_fun*vv;
   inv_Jts = inv(Jts);
   det_Jts = abs(det(Jts));

   Gs = inv_Jts*interp_fun;

   %int = 0;
   %for q=1:4
   %  int = int + w(q)*Gs'*Gs;
   %end
   %Inte = int*det_Jts; % or simplified
   % Citation: The FEM displayed G.Dhatt & G. Touzot

   Inte = (1/6)*Gs'*Gs*det_Jts;

   Local = sparse(vr,vr);
   Local(simp(i,:),simp(i,:))= Inte;
   IntGrad(:,i) = Local(:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
