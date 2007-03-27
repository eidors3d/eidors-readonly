function [srf] = find_boundary(simp);
% [srf] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  Outter boundary surfaces
%simp = The simplices matrix

wew = size(simp,2) - 1;

S=[];

for b=1:size(simp,1)
   x = simp(b,1);
   y = simp(b,2);
   z = simp(b,3);
   w = simp(b,4);
   s1 = [x y z];
   s2 = [x y w];
   s3 = [y z w];
   s4 = [x z w];
   
   Sn = [s1; s2; s3; s4];
   S =  [S ; Sn];
end

N=sort(S,2);
M=sortrows(N);

count=1;
inc=1;

m=size(M,1);
for i=1:count:m-1
   
   ithrow = M(i,:);
   jthrow = M(i+1,:); 
   
   if ithrow == jthrow 
      M(i,:)=zeros(1,wew);
      M(i+1,:)=zeros(1,wew);
      count=2;
   end
   
   if ithrow ~= jthrow
   count=1;
   end
   
 end

for dk=1:m

if M(dk,:)~= zeros(1,wew);
     XL(inc,:)= M(dk,:);
      inc=inc+1;
  end
   
end
XL;
srf = XL;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
