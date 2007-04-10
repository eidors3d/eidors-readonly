function [srf] = find_boundary(simp);
% [srf] = find_boundary(simp);
%
%Caclulates the boundary faces of a given 3D volume.
%Usefull in electrode assignment.
%
%srf  =  array of elements on each boundary simplex
%        boundary simplices are of 1 lower dimention than simp
%simp = The simplices matrix

% $Id: find_boundary.m,v 1.3 2007-04-10 14:40:28 aadler Exp $

wew = size(simp,2) - 1;

if wew==3 || wew==2
   srf= find_2or3d_boundary(simp,wew);
else
   error('not 2D or 3D simplices');
end

function srf= find_2or3d_boundary(simp,wew);
els= size(simp,1);
srf= [];
% find elements on the boundary
ks= zeros(els,wew+1);
for i=1:els
   for j=1:wew+1
      ks(:,j)= sum(simp== simp(i,j),2);
   end
   fs= (  sum(ks,2) == wew  );
   if sum(fs) < wew+1; % elem is on boundary 
      % sel is number of times each point is seen
      sel = sum(ks(fs,:),1);
      % if elem is on one boundary  , sel has one element  == wew
      % if elem is on two boundaries, sel has two elements == wew-1
      m_sel = max(sel);
      for ff= find( sel == m_sel );
         idx = 1:wew+1;
         idx(ff) = [];
         srf= [srf; simp(i,idx)]; 
      end
   end
end

% sort the output srf
   srf = sort(srf,2);
   srf = sortrows(srf);

function srf= findk3d_boundary(simp,wew);
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
