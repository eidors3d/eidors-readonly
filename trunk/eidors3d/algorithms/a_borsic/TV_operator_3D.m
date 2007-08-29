function L=TV_operator_3D( msh )
% L=TV_operator_3D( msh )
% construct Total Variation operator for a 3D mesh
%
% INPUTS:
% msh = a 3D scaip msh structure with .msh.vtx_c_c, .elem_c defined
% OUTPUTS:
% L = TV operator (generally used for regularisation)
%
% Copyright 2004 Andrea Borsic, SC-AIP s.a.s.
% Scientific Computing & Applied Inverse Problems  www.sc-aip.com
% Modifications (C) 2005 Andy Adler.
% License: GPL version 2 or version 3
% $Id: TV_operator_3D.m,v 1.4 2007-08-29 09:13:50 aadler Exp $

%global SP3D % Sparse 3D matrix used in the computations
%SP3D=[];

num_vtx=size(msh.vtx_c,1);
num_tet=size(msh.elem_c,1);

list=[]; % List is an auxiliary variable which will hold for each row the two facing thetrahydra and the shared face area

%We build a selection array, to index the T matrix

SEL=[ 1 1 1 2; 2 3 4 4; 3 4 2 3];

% allocate_3Dsparse(num_vtx,num_vtx,num_vtx,num_tet*4);
SP3D=spalloc(num_vtx,num_vtx^2,num_tet*4);

for k=1:num_tet
    
    for j=1:4 % cycle on the faces of each thetrahydra
        
        face_vtx(1)=msh.elem_c(k,SEL(1,j)); % face_vtx are the varticies on a face of thetrahydra k
        face_vtx(2)=msh.elem_c(k,SEL(2,j));
        face_vtx(3)=msh.elem_c(k,SEL(3,j));
        
        face_vtx=sort(face_vtx); % faces must be unique
        
%       simplex=read_from_3Dsparse(face_vtx,num_vtx,num_vtx,num_vtx);
%function val=read_from_3Dsparse(mnp,M,N,P);
%
%global SP3D
        %val=SP3D(mnp(1),(mnp(2)-1)*N+mnp(3));
        simplex=SP3D(face_vtx(1),(face_vtx(2)-1)*num_vtx+face_vtx(3));
        
        if (simplex==0)
            
%         write_to_3Dsparse(mnp,M,N,P,val);
%           SP3D(mnp(1),(mnp(2)-1)*N+mnp(3))=val;
%
%           write_to_3Dsparse(face_vtx,num_vtx,num_vtx,num_vtx,k);
            SP3D(face_vtx(1),(face_vtx(2)-1)*num_vtx+face_vtx(3))=k;
            
        else
            
            vec1=msh.vtx_c(face_vtx(1),:)-msh.vtx_c(face_vtx(2),:);  % vec1,vec2 are vectors along edgs, used for the area calculation as cross prod.
            vec2=msh.vtx_c(face_vtx(1),:)-msh.vtx_c(face_vtx(3),:);
            facearea=0.5*norm(cross(vec1,vec2));
            list=[list;[k,simplex,facearea]]; % Triangles and length of the shared edge are written into the list
            
        end % if then else
        
    end % for j
    
end % for k

L=spalloc(length(list),num_tet,2*length(list));

for i=1:length(list)
    
    L(i,list(i,1))=list(i,3);
    L(i,list(i,2))=-list(i,3);
    
end % for

%clear SP3D;


%%%%%%%%%% Auxiliary functions for handling the 3D sparse matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% mnp is a vector of 3 indexes into the matrix, and M,N,P the matrix size along the three dimensions
% - removed these functions for inclusion with EIDORS -aa Dec05
%
%function allocate_3Dsparse(M,N,P,max_elements);
%
%global SP3D
%SP3D=spalloc(M,N*P,max_elements);
%
%function write_to_3Dsparse(mnp,M,N,P,val);
%
%global SP3D
%SP3D(mnp(1),(mnp(2)-1)*N+mnp(3))=val;
%
%function val=read_from_3Dsparse(mnp,M,N,P);
%
%global SP3D
%val=SP3D(mnp(1),(mnp(2)-1)*N+mnp(3));
%
