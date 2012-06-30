function[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(filename)
% Function to read in a mesh model from Gmsh and saves it in
% five arrays; surface (srf), veritices (vtx), face no. (fc)
% volume (simp) and edges (edg)
%
% srf        = The surfaces indices into vtx
% simp       = The volume indices into vtx
% vtx        = The vertices matrix
% fc         = A one column matrix containing the face numbers
% edg        = Edge segment information
% filename   = Name of file containing NetGen .vol information
% mat_ind    = Material index
% phys_names = Structure of "Physical" entities in the mesh
%              .dim   = dimension
%              .name  = name (string)
%              .tag   = physical tag
%              .nodes = N-x-dim array of indices into vtx

% $Id$
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2
% Modified by James Snyder <jbsnyder@fanplastic.org>

fid = fopen(filename,'r');
phys_names = [];
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); break; end

    if strcmp(tline,'$Elements')
       elements= parse_elements( fid );
    elseif strcmp(tline,'$Nodes')
       nodes= get_lines_with_nodes( fid );
    elseif strcmp(tline,'$PhysicalNames')
       phys_names= parse_names( fid );
    end
end

if ~isempty(phys_names)
    for i = 1:length(phys_names)
        tmpelements = find(arrayfun(@(x)x.phys_tag==phys_names(i).tag,elements));
        phys_names(i).nodes = cat(1,elements(tmpelements).simp);
    end
end

edg = [];
bc = [];
mat_ind = [];

% Select 2d vs 3d model by checking if Z is all the same
if length( unique( nodes(:,4) ) ) > 1 
    vtx = nodes(:,2:4);
    % Type 2: 3-node triangle
    tri = find(arrayfun(@(x)x.type==2,elements));
    % Type 4: 4-node tetrahedron
    tet = find(arrayfun(@(x)x.type==4,elements));
    simp = cat(1,elements(tet).simp);
    srf = cat(1,elements(tri).simp);
else
    vtx = nodes(:,2:3);
    tri = find(arrayfun(@(x)x.type==2,elements));
    simp = cat(1,elements(tri).simp);
    srf = [];
end

elemtags = cat(1,elements.phys_tag);
fc = elemtags(tri,1);
end


function mat = get_lines_with_nodes( fid )
% Line Format:
% node-number x-coord y-coord z-coord
tline = fgetl(fid);
n_rows = sscanf(tline,'%d');
mat= fscanf(fid,'%f',[4,n_rows])';
end

function names = parse_names( fid )
% Line Format:
% physical-dimension physical-number "physical-name"
tline = fgetl(fid);
n_rows = sscanf(tline,'%d');
names = struct('tag',{},'dim',{},'name',{});
for i = 1:n_rows
    tline = fgetl(fid);
    if exist('OCTAVE_VERSION')
        parts = strsplit(tline,' ');
    else
        parts = regexp(tline,' ','split');
    end
    nsz = size(names,2)+1;
    names(nsz).dim = str2double( parts(1) );
    names(nsz).tag = str2double( parts(2) );
    tname = parts(3);
    names(nsz).name = strrep(tname{1},'"','');
end
end

function elements = parse_elements( fid )
% Line Format:
% elm-number elm-type number-of-tags < tag > ... node-number-list
tline = fgetl(fid);
n_rows = sscanf(tline,'%d');
elements = struct('simp',{},'phys_tag',{},'geom_tag',{});
for i = 1:n_rows
    tline = fgetl(fid);
    n = sscanf(tline, '%d')';
    nsz = size(elements,2)+1;
    elements(nsz).simp = n(n(3) + 4:end);
    % 
    elements(nsz).type = n(2);
    if n(3) > 0 % get tags if they exist
        % By default, first tag is number of parent physical entity
        % second is parent elementary geometrical entity
        % third is number of parent mesh partitions followed by
        % partition ids
        tags = n(4:3+n(3));
        if length(tags) >= 1
            elements(nsz).phys_tag = tags(1);
            if length(tags) >= 2
                elements(nsz).geom_tag = tags(2);
            end
        end
    end
end
end
