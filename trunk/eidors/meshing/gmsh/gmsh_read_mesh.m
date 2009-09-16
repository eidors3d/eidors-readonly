function[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(filename)
% Function to read in a mesh model from NetGen and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
% 
% Version 4.0
% B.D.Grieve - 27/01/2002 + modifications by lmazurk
% A.Adler - 2006 mods to run quicker
%
% EIDORS's srf array is a subset of NetGen's surface element data
% (columns 6:8). The first column of the surface element data also
% ascribes a face number to each surface which is saved as the fc
% array. Each line of the srf array contains 3 indices to define
% a triangle mapped on to the three dimensional vtx array.
% EIDORS's vtx array is a direct equivalent to NetGen's pointer data.
%
%
% srf      = The surfaces indices into vtx
% simp     = The volume indices into vtx
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% edg      = Edge segment information
% filename = Name of file containing NetGen .vol information
% mat_ind  = Material index

% $Id: ng_read_mesh.m 1535 2008-07-26 15:36:27Z aadler $
% (C) 2002-2006 (C) Licenced under the GPL


fid = fopen(filename,'r');
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); break; end

    if strcmp(tline,'$Elements')
       elements= get_lines_with_elements( fid );
    elseif strcmp(tline,'$Nodes')
       nodes= get_lines_with_nodes( fid );
    end
end

%srf = elements(:,6:8);
srf = [];
%fc = elements(:,1);
fc = [];
simp = elements(:,7:9);
%edg = elements;
edg = [];
mat_ind=elements(:,5);
%bc = elements(:,2);
bc = [];
vtx = nodes(:,2:3);
end

function mat= get_lines_with_nodes( fid )
   tline = fgetl(fid);
   n_rows = sscanf(tline,'%d');
   mat= fscanf(fid,'%f',[4,n_rows])';
end

function mat = get_lines_with_elements( fid )
   tline = fgetl(fid);
   n_rows = sscanf(tline,'%d');
   mat = [];
   for i = 1:n_rows
       tline = fgetl(fid);
       parts = regexp(tline,' ','split');
       if str2double( parts(2) ) == 2
           n = str2double( parts );
           mat(size(mat,1)+1,:) = n;
       end
   end
end
