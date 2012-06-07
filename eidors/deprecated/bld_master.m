function [Ef,D,Ela] = bld_master(vtx,simp,mat_ref);
%function [Ef,D,Ela] = bld_master(vtx,simp,mat_ref);
%
%Builds up the main compartment (GAP-SHUNT) of the system matrix 
%for the complete electrode model. It is called within the function 
%fem_master_full.
%
%
%
%Ef      = The UNreferenced GAP based system matrix
%D       = The sgradients of the shape functions 
%          over each element.
%Ela     = Normalised volums of the elements
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%mat_ref = The reference CONDUCTIVITY at each element. 
%In the complex case mat_ref(i) = sigma(i) - epsilon(i)

warning('EIDORS:deprecated','BLD_MASTER is deprecated as of 07-Jun-2012. ');

dimen= size(vtx,2);
if dimen==2
   [Ef,D,Ela] = bld_master_2d(vtx,simp,mat_ref); 
elseif dimen==3
   [Ef,D,Ela] = bld_master_3d(vtx,simp,mat_ref);
else
   error('not 2d or 3d');
end
   
function [Ef,A,Ela] = bld_master_2d(vtx,simp,mat_ref) 
%Based on gm_assemble by S.Vavasis , http://www.cs.cornel.edu/ 

[vr, vc] = size(vtx);
[sr, sc] = size(simp);
a = mat_ref;

ilist = kron((1:vc*sr), [1,1]);
jlist = zeros(1,sr*vc*2);
slist = zeros(1,sr*vc*2);

for d = 1 : vc
  jlist(2 * d - 1 : 2 * vc : sr * vc * 2) = simp(:,1)';
  jlist(2 * d : 2 * vc : sr * vc * 2) = simp(:, d + 1);
  slist(2 * d - 1 : 2 * vc : sr * vc * 2) = -ones(1,sr);
  slist(2 * d : 2 * vc : sr * vc * 2) = ones(1,sr);
end

A0 = sparse(ilist,jlist,slist,vc*sr,vr);


if vc == 2
  J1 = A0 * vtx(:,1);
  J2 = A0 * vtx(:,2);
  J11 = J1(1:2:sr*2);
  J12 = J1(2:2:sr*2);
  J21 = J2(1:2:sr*2);
  J22 = J2(2:2:sr*2);
  detJ = J11 .* J22 - J21 .* J12;
  
  
  invJ11 = J22 ./ detJ;
  invJ12 = -J12 ./ detJ;
  invJ21 = -J21 ./ detJ;
  invJ22 = J11 ./ detJ;
elseif vc == 3
  J1 = A0 * vtx(:,1);
  J2 = A0 * vtx(:,2);
  J3 = A0 * vtx(:,3);
  J11 = J1(1:3:sr*3);
  J12 = J1(2:3:sr*3);
  J13 = J1(3:3:sr*3);
  J21 = J2(1:3:sr*3);
  J22 = J2(2:3:sr*3);
  J23 = J2(3:3:sr*3);
  J31 = J3(1:3:sr*3);
  J32 = J3(2:3:sr*3);
  J33 = J3(3:3:sr*3);
  detJ = J11 .* J22 .* J33 - J11 .* J23 .* J32 - J12 .* J21 .* J33 ...
          + J12 .* J23 .* J31 + J13 .* J21 .* J32 - J13 .* J22 .* J31;
       
       
  invJ11 = (J22 .* J33 - J23 .* J32) ./ detJ;
  invJ12 = (J32 .* J13 - J12 .* J33) ./ detJ;
  invJ13 = (J12 .* J23 - J22 .* J13) ./ detJ;
  invJ21 = (J31 .* J23 - J21 .* J33) ./ detJ;
  invJ22 = (J11 .* J33 - J13 .* J31) ./ detJ;
  invJ23 = (J21 .* J13 - J11 .* J23) ./ detJ;
  invJ31 = (J21 .* J32 - J31 .* J22) ./ detJ;
  invJ32 = (J31 .* J12 - J11 .* J32) ./ detJ;
  invJ33 = (J11 .* J22 - J21 .* J12) ./ detJ;
else
  error('Master matrix construction failed')
end


ilist = kron((1 : vc * sr), ones(1,vc));
jlist = zeros(1, sr*vc^2);
for d = 1 : vc 
  jlist(d:vc:sr*vc^2) = kron((d:vc:vc*sr),ones(1,vc));
end

if vc == 2
  slist = zeros(1,sr*4);
  slist(1:4:sr*4) = invJ11;
  slist(2:4:sr*4) = invJ21;
  slist(3:4:sr*4) = invJ12;
  slist(4:4:sr*4) = invJ22;
else
  slist = zeros(1,sr*9);
  slist(1:9:sr*9) = invJ11;
  slist(2:9:sr*9) = invJ21;
  slist(3:9:sr*9) = invJ31;
  slist(4:9:sr*9) = invJ12;
  slist(5:9:sr*9) = invJ22;
  slist(6:9:sr*9) = invJ32;
  slist(7:9:sr*9) = invJ13;
  slist(8:9:sr*9) = invJ23;
  slist(9:9:sr*9) = invJ33;
end


ElJac = sparse(ilist,jlist,slist,vc*sr,vc*sr);
A = ElJac * A0;


Vols = abs(detJ) / prod(1:vc);

materials = length(a);
volumes = size(Vols);

if materials ~= volumes
  error('Some elements have not been assigned');
end;

Ela = sparse( (1:vc*sr), (1:vc*sr),kron( (a .* Vols)',ones(1,vc)) );

Ef = A'*Ela*A; 

%This is for the Jacobian matrix (does not include conductivity)
Ela = sparse( (1:vc*sr), (1:vc*sr),kron( Vols.',ones(1,vc)) );


%
% 3D BLD MASTER
%
function [Ef,D,Ela] = bld_master_3d(vtx,simp,mat_ref) 
[nv, dimen] = size(vtx);  %Number of vertices and dimension
[ns, dimen_p1] = size(simp); %Number of simplices
a = mat_ref;
dimen2 = 2*dimen;

ils = 1:dimen*ns;
ilst(1:2:dimen2*ns) = ils;
ilst(2:2:dimen2*ns) = ils;


patt = 1:dimen2:ns*dimen2;

jlst(patt) = simp(:,1);
jlst(patt+1) = simp(:,2);
jlst(patt+2) = simp(:,1);
jlst(patt+3) = simp(:,3);
jlst(patt+4) = simp(:,1);
jlst(patt+5) = simp(:,4);


sls = ones(1,dimen*ns);
slst(1:2:dimen2*ns) = -sls;
slst(2:2:dimen2*ns) = sls;


D0 = sparse(ilst,jlst,slst,dimen*ns,nv); 
%D0 is the matrix of the definitions of the gradients on elements

J1 = D0 * vtx(:,1);
J2 = D0 * vtx(:,2);
J3 = D0 * vtx(:,3);
  
JJ= zeros( 3,3, ns );
for w=1:3
   r=1;
   JJ(r,w,:)   = J1(w:dimen:ns*dimen);
   JJ(r+1,w,:) = J2(w:dimen:ns*dimen);
   JJ(r+2,w,:) = J3(w:dimen:ns*dimen);
end

dj = squeeze(sum( [prod([JJ(1,1,:);JJ(2,2,:);JJ(3,3,:)]); prod([JJ(1,2,:);JJ(2,3,:);JJ(3,1,:)]);...
                   prod([JJ(1,3,:);JJ(2,1,:);JJ(3,2,:)]); prod([-JJ(1,3,:);JJ(2,2,:);JJ(3,1,:)]);...
                   prod([-JJ(1,1,:);JJ(2,3,:);JJ(3,2,:)]); prod([-JJ(1,2,:);JJ(2,1,:);JJ(3,3,:)])]));

      
ilst = kron((1:dimen*ns), ones(1,dimen));
jlst = zeros(1, ns*dimen^2);
for d = 1:dimen 
  jlst(d:dimen:ns*dimen^2) = kron((d:dimen:dimen*ns),ones(1,dimen));
end
slst = zeros(1,ns*dimen^2); 
   
  pat = 1:dimen^2:ns*dimen^2; 
  
  slst(pat) = squeeze(sum([prod([JJ(2,2,:) ; JJ(3,3,:)]); prod([-JJ(2,3,:) ; JJ(3,2,:)])])) ./ dj; 
  slst(pat+1) = squeeze(sum([prod([JJ(3,1,:) ; JJ(2,3,:)]); prod([-JJ(2,1,:) ; JJ(3,3,:)])])) ./ dj; 
  slst(pat+2) = squeeze(sum([prod([JJ(2,1,:) ; JJ(3,2,:)]); prod([-JJ(3,1,:) ; JJ(2,2,:)])])) ./ dj; 
  slst(pat+3) = squeeze(sum([prod([JJ(3,2,:) ; JJ(1,3,:)]); prod([-JJ(1,2,:) ; JJ(3,3,:)])])) ./ dj; 
  slst(pat+4) = squeeze(sum([prod([JJ(1,1,:) ; JJ(3,3,:)]); prod([-JJ(1,3,:) ; JJ(3,1,:)])])) ./ dj; 
  slst(pat+5) = squeeze(sum([prod([JJ(3,1,:) ; JJ(1,2,:)]); prod([-JJ(1,1,:) ; JJ(3,2,:)])])) ./ dj; 
  slst(pat+6) = squeeze(sum([prod([JJ(1,2,:) ; JJ(2,3,:)]); prod([-JJ(2,2,:) ; JJ(1,3,:)])])) ./ dj; 
  slst(pat+7) = squeeze(sum([prod([JJ(2,1,:) ; JJ(1,3,:)]); prod([-JJ(1,1,:) ; JJ(2,3,:)])])) ./ dj; 
  slst(pat+8) = squeeze(sum([prod([JJ(1,1,:) ; JJ(2,2,:)]); prod([-JJ(2,1,:) ; JJ(1,2,:)])])) ./ dj; 
  

LocJac = sparse(ilst,jlst,slst,dimen*ns,dimen*ns);

D = LocJac * D0;


Vols = abs(dj(:)) / 6;

materials = length(a);
volumes = size(Vols);

if materials ~= volumes
  error('Some elements have not been assigned');
end;

%This is for the global conductance matrix (includes conductivity)
Ela = sparse( (1:dimen*ns), (1:dimen*ns),kron( (a.* Vols).',ones(1,dimen)) );

Ef = D'*Ela*D; 

%This is for the Jacobian matrix (does not include conductivity)
Ela = sparse( (1:dimen*ns), (1:dimen*ns),kron( Vols.',ones(1,dimen)) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of the EIDORS suite.
% Copyright (c) N. Polydorides 2003
% Copying permitted under terms of GNU GPL
% See enclosed file gpl.html for details.
% EIDORS 3D version 2.0
% MATLAB version 5.3 R11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

