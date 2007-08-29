function [D,Vols] = calc_shape_grads(fwd_model);
%function [Y,D,Y_homo,Vols] = bld_master_PointElectrodes(vtx,simp,mat_ref);
%Builds up the main compartment (GAP-SHUNT) of the system matrix 
%for the complete electrode model. It is called within the function 
%fem_master_full.
%
%D       = The gradients of the shape functions over each element.
%Vols     = Normalised volumes of the elements
%
% (C) 2005  Brad Graham.  %  Licenced under GNU GPL
% $Id: calc_shape_grads.m,v 1.9 2007-08-29 09:20:55 aadler Exp $

NODE = fwd_model.nodes;
ELEM = fwd_model.elems;
    
[nv, dimen] = size(NODE);  %Number of vertices and dimension
[ns, dimen_p1] = size(ELEM); %Number of simplices (elements)
dimen2 = 2*dimen;

ils = 1:dimen*ns;
ilst(1:2:dimen2*ns) = ils;
ilst(2:2:dimen2*ns) = ils;

patt = 1:dimen2:ns*dimen2;

jlst(patt) = ELEM(:,1);
jlst(patt+1) = ELEM(:,2);
jlst(patt+2) = ELEM(:,1);
jlst(patt+3) = ELEM(:,3);

if(dimen==3)
    % disp('3D.');
    jlst(patt+4) = ELEM(:,1);  %3D
    jlst(patt+5) = ELEM(:,4);  %3D
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
    a2=NODE(ELEM(ii,:),:);
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
    ii=ELEM(jj,:);
    if(dimen==3)
        ii=[ii(2:4) ii(1)]; %3D
    else
        ii=[ii(2:3) ii(1)]; % 2D    
    end
    nodes_ = NODE(ii,:);
    J_k=M_0*nodes_;
    invJ_k=inv(J_k)';
    
    tw=invJ_k(:)';
    slst=[slst, tw];
end

%LocJac = sparse(ilst,jlst,slst,dimen*ns,dimen*ns);
LocJac = sparse(ilst,jlst,slst,dimen*ns,dimen*ns).*sqrt(2);
%Brad's correction to make this code match Silvestri's

D = LocJac * D0;

Vols = abs(dj) / factorial( dimen);


