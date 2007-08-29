function [Y,D,Y_homo,Vols] = bg_bld_master(mesh,param2, param3);
%function [Y,D,Y_homo,Vols] = bld_master_PointElectrodes(vtx,simp,mat_ref);
%Builds up the main compartment (GAP-SHUNT) of the system matrix 
%for the complete electrode model. It is called within the function 
%fem_master_full.
%
%Y      = The UNreferenced GAP based system matrix
%D       = The gradients of the shape functions 
%          over each element.
%Vols     = Normalised volumes of the elements
%vtx     = The vertices matrix.
%simp    = The simplices matrix.
%mat_ref = The reference CONDUCTIVITY at each element. 
%In the complex case mat_ref(i) = sigma(i) - epsilon(i)
%
% (C) 2003 Nick Polydorides, modified by  Brad Graham (2005)
%  Licenced under GNU GPL
% $Id: bg_bld_master.m,v 1.4 2007-08-29 09:13:50 aadler Exp $

if nargin == 3
    mesh.NODE = mesh;
    mesh.ELEM = param2;
    mesh.sigma = param3;
elseif nargin ~=1
    error('wrong number of args');
end
    
[nv, dimen] = size(mesh.NODE);  %Number of vertices and dimension
[ns, dimen_p1] = size(mesh.ELEM); %Number of simplices (elements)
dimen2 = 2*dimen;

ils = 1:dimen*ns;
ilst(1:2:dimen2*ns) = ils;
ilst(2:2:dimen2*ns) = ils;

patt = 1:dimen2:ns*dimen2;

jlst(patt) = mesh.ELEM(:,1);
jlst(patt+1) = mesh.ELEM(:,2);
jlst(patt+2) = mesh.ELEM(:,1);
jlst(patt+3) = mesh.ELEM(:,3);

if(dimen==3)
    % disp('3D.');
    jlst(patt+4) = mesh.ELEM(:,1);  %3D
    jlst(patt+5) = mesh.ELEM(:,4);  %3D
else
    % disp('2D.');
end

sls = ones(1,dimen*ns);
slst(1:2:dimen2*ns) = -sls;
slst(2:2:dimen2*ns) = sls;

D0 = sparse(ilst,jlst,slst,dimen*ns,nv); 
%D0 is the matrix of the definitions of the gradients on elements


if(dimen==3)
    aa=[1 1 1 1]'; %3D
else
    aa=[1 1 1]';    %2D
end

for ii=1:ns
    a2=mesh.NODE(mesh.ELEM(ii,:),:);
    a2=[aa a2];
    dj2(ii)=det(a2);
end
dj=dj2';

ilst = kron((1:dimen*ns), ones(1,dimen));

jlst = zeros(1, ns*dimen^2);
for d = 1:dimen 
    jlst(d:dimen:ns*dimen^2) = kron((d:dimen:dimen*ns),ones(1,dimen));
end


if(dimen==3)
    M_0 = [1 0 0 -1;0 1 0 -1; 0 0 1 -1]; % 3D
else
    M_0 = [1 0 -1;0 1 -1];    %2D    
end

slst=[];
for jj=1:ns
    ii=mesh.ELEM(jj,:);
    if(dimen==3)
        ii=[ii(2:4) ii(1)]; %3D
    else
        ii=[ii(2:3) ii(1)]; % 2D    
    end
    nodes_ = mesh.NODE(ii,:);
    J_k=M_0*nodes_;
    invJ_k=inv(J_k)';
    
    tw=invJ_k(:)';
    slst=[slst, tw];
end

%LocJac = sparse(ilst,jlst,slst,dimen*ns,dimen*ns);
LocJac = sparse(ilst,jlst,slst,dimen*ns,dimen*ns).*sqrt(2);
%Brad's correction to make this code match Silvestri's

D = LocJac * D0;

if (dimen==3)
    Vols = abs(dj) / 6; %3D
else %dimen==2
    Vols = abs(dj) / 2; %2D
end

if length(mesh.sigma) ~= size(Vols)
    error('Some elements have not been assigned');
end;

%This is for the global conductance matrix (includes conductivity)
Ela = sparse( (1:dimen*ns), (1:dimen*ns),kron( (mesh.sigma.* Vols).',ones(1,dimen)) );

Y = D'*Ela*D; 

%This is for the Jacobian calculation function, matrix (does not include conductivity)
Y_homo = sparse( (1:dimen*ns), (1:dimen*ns),kron( Vols.',ones(1,dimen)) );

