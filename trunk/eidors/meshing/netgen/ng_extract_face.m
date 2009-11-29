function[fcsrf,fci] = ng_extract_face(srf,vtx,fc,fcnmb)

% This function takes the wireframe model created by readmesh
% and extracts one of the faces (fcnmb) into an indices
% matrix (fci) mapped on to indices matrix srf which, in turn,
% maps on to vtx. It also creates a matrix (fcsrf) which
% directly maps this face into vtx.
%
% Version 2.4
% B.D.Grieve - 23/01/2002
%
%
% srf      = The boundary surfaces
% vtx      = The vertices matrix
% fc       = The face numbers to which each surface belong
% fcnmb    = The number of the face to be extracted
% fcsrf    = The indices into the vtx matrix for this face
% fci      = The indices into the srf matrix of this face

 fci  = find( fc == fcnmb );
 fcsrf= srf(fci,:);

 return 
% OLDER, COMPLICATED CODE
fcsrf = [];
fci = [];
for loop1 = 1:size(fc,1)
    if fc(loop1)==fcnmb
        fcsrf = [fcsrf; srf(loop1,:)]; % Concatenate fcsrf to map directly with vtx for this face
        fci = [fci; loop1]; % Concatenate fci with srf indices for this face
    end
end
