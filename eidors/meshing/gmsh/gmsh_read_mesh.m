function[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh(filename)
% Function to read in a mesh model from Gmsh and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
%
% srf      = The surfaces indices into vtx
% simp     = The volume indices into vtx
% vtx      = The vertices matrix
% fc       = A one column matrix containing the face numbers
% edg      = Edge segment information
% filename = Name of file containing NetGen .vol information
% mat_ind  = Material index

% $Id$
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2


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

surf_mesh = elements(:,end) ==0;
if any(surf_mesh)
   srf = elements(surf_mesh,6:end-1);
   elements(surf_mesh) = [];
end

%fc = elements(:,1);
fc = [];
simp = elements(:,6:end);
%edg = elements;
edg = [];
mat_ind=elements(:,5);
%bc = elements(:,2);
bc = [];
vtx = nodes(:,2:end);
end

function mat= get_lines_with_nodes( fid )
   tline = fgetl(fid);
   n_rows = sscanf(tline,'%d');
   mat= fscanf(fid,'%f',[4,n_rows])';
end

function mat = get_lines_with_elements( fid )
   tline = fgetl(fid);
   n_rows = sscanf(tline,'%d');
   mat = zeros(n_rows,8);
   count = 1;
   for i = 1:n_rows
       tline = fgetl(fid);
       parts = sscanf(tline,'%f');
       if parts(2) == 2 || parts(2) == 4
           mat(count,1:length(parts)) = parts;
           count = count + 1;
       end
   end
   if count <= n_rows
      mat(count:end,:) = []; % remove unused rows
   end
end
