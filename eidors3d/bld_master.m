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

