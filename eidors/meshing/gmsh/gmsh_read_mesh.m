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

if ischar(filename) && strcmp(filename,'UNIT_TEST'); do_unit_test; return; end

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
% elements = struct('simp',{},'phys_tag',{},'geom_tag',{});
elements(n_rows).simp = [];
elements(n_rows).phys_tag = [];
elements(n_rows).geom_tag = [];
elements(n_rows).type = [];

for i = 1:n_rows
    tline = fgetl(fid);
    n = sscanf(tline, '%d')';
%     nsz = size(elements,2)+1;
    elements(i).simp = n(n(3) + 4:end);
    % 
    elements(i).type = n(2);
    if n(3) > 0 % get tags if they exist
        % By default, first tag is number of parent physical entity
        % second is parent elementary geometrical entity
        % third is number of parent mesh partitions followed by
        % partition ids
        tags = n(4:3+n(3));
        if length(tags) >= 1
            elements(i).phys_tag = tags(1);
            if length(tags) >= 2
                elements(i).geom_tag = tags(2);
            end
        end
    end
end
end

function do_unit_test
   tmpnam = tempname;
   fid = fopen(tmpnam,'w');
   fprintf(fid,gmshv2file);
   fclose(fid);
   [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(tmpnam);
    unit_test_cmp('v2 vtx ',vtx(2:3,:),[1,0;-1,0])
    unit_test_cmp('v2 simp',simp(2:3,:),[2,4,15; 14,17,19]);

   fid = fopen(tmpnam,'w');
   fprintf(fid,gmshv4file);
   fclose(fid);
   [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(tmpnam);
    
end

% Example of gmsh v4 file
function t = gmshv4file; t=[ ...
'$MeshFormat\n' ...
'4 0 8\n' ...
'$EndMeshFormat\n' ...
'$Entities\n' ...
'3 2 1 0\n' ...
'1 0 0 0 0 0 0 0 \n' ...
'2 1 0 0 1 0 0 0 \n' ...
'3 -1 0 0 -1 0 0 0 \n' ...
'1 -1 0 0 1 0.984807753012208 0 0 2 2 -3 \n' ...
'2 -1 -0.984807753012208 0 1 0 0 0 2 3 -2 \n' ...
'4 -1 -0.984807753012208 0 1 0.984807753012208 0 0 2 1 2 \n' ...
'$EndEntities\n' ...
'$Nodes\n' ...
'6 14\n' ...
'1 0 0 1\n' ...
'1 0 0 0\n' ...
'2 0 0 1\n' ...
'2 1 0 0\n' ...
'3 0 0 1\n' ...
'3 -1 0 0\n' ...
'1 1 0 3\n' ...
'4 0.7071067795891069 0.7071067827839882 0\n' ...
'5 -4.62506049552981e-09 1 0\n' ...
'6 -0.7071067828954029 0.7071067794776923 0\n' ...
'2 1 0 3\n' ...
'7 -0.7071067795891069 -0.7071067827839882 0\n' ...
'8 4.62506049552981e-09 -1 0\n' ...
'9 0.7071067828954029 -0.7071067794776923 0\n' ...
'4 2 0 5\n' ...
'10 -7.225004064876332e-18 5.806922560102616e-17 0\n' ...
'11 0.4208610471654206 -0.174326352879787 0\n' ...
'12 -0.4208610471654206 0.1743263528797871 0\n' ...
'13 0.1601886191568594 0.3867295406713929 0\n' ...
'14 -0.1601886191568594 -0.3867295406713928 0\n' ...
'$EndNodes\n' ...
'$Elements\n' ...
'6 27\n' ...
'1 0 15 1\n' ...
'1 1 \n' ...
'2 0 15 1\n' ...
'2 2 \n' ...
'3 0 15 1\n' ...
'3 3 \n' ...
'1 1 1 4\n' ...
'4 2 4 \n' ...
'5 4 5 \n' ...
'6 5 6 \n' ...
'7 6 3 \n' ...
'2 1 1 4\n' ...
'8 3 7 \n' ...
'9 7 8 \n' ...
'10 8 9 \n' ...
'11 9 2 \n' ...
'4 2 2 16\n' ...
'12 2 4 11 \n' ...
'13 3 7 12 \n' ...
'14 12 7 14 \n' ...
'15 11 4 13 \n' ...
'16 5 6 12 \n' ...
'17 8 9 11 \n' ...
'18 8 11 14 \n' ...
'19 5 12 13 \n' ...
'20 10 11 13 \n' ...
'21 10 12 14 \n' ...
'22 4 5 13 \n' ...
'23 7 8 14 \n' ...
'24 12 10 13 \n' ...
'25 11 10 14 \n' ...
'26 6 3 12 \n' ...
'27 9 2 11 \n' ...
'$EndElements\n'];
end

% Example of gmsh v2 file
function t = gmshv2file; t=[ ...
'$MeshFormat\n' ...
'2.2 0 8\n' ...
'$EndMeshFormat\n' ...
'$Nodes\n' ...
'26\n' ...
'1 0 0 0\n' ...
'2 1 0 0\n' ...
'3 -1 0 0\n' ...
'4 0.8660254037855228 0.4999999999981222 0\n' ...
'5 0.500000000002617 0.8660254037829277 0\n' ...
'6 3.95800631512864e-12 1 0\n' ...
'7 -0.4999999999980611 0.866025403785558 0\n' ...
'8 -0.8660254037843073 0.5000000000002276 0\n' ...
'9 -0.8660254037855228 -0.4999999999981222 0\n' ...
'10 -0.500000000002617 -0.8660254037829277 0\n' ...
'11 -3.95800631512864e-12 -1 0\n' ...
'12 0.4999999999980611 -0.866025403785558 0\n' ...
'13 0.8660254037843073 -0.5000000000002276 0\n' ...
'14 5.368222989669768e-17 1.262108679519644e-17 0\n' ...
'15 0.5000000000000001 0.1339745962149803 0\n' ...
'16 -0.5 -0.1339745962149803 0\n' ...
'17 0.133974596217502 0.5000000000002625 0\n' ...
'18 -0.133974596217502 -0.5000000000002626 0\n' ...
'19 -0.3660254037841135 0.3660254037850291 0\n' ...
'20 0.3660254037841134 -0.3660254037850291 0\n' ...
'21 0.6830127018922291 -0.1830127018925211 0\n' ...
'22 -0.683012701892229 0.1830127018925211 0\n' ...
'23 0.5000000000014104 0.4999999999990732 0\n' ...
'24 -0.5000000000014104 -0.4999999999990732 0\n' ...
'25 -0.1830127018901805 0.6830127018929458 0\n' ...
'26 0.1830127018901805 -0.6830127018929458 0\n' ...
'$EndNodes\n' ...
'$Elements\n' ...
'51\n' ...
'1 15 2 0 1 1\n' ...
'2 15 2 0 2 2\n' ...
'3 15 2 0 3 3\n' ...
'4 1 2 0 1 2 4\n' ...
'5 1 2 0 1 4 5\n' ...
'6 1 2 0 1 5 6\n' ...
'7 1 2 0 1 6 7\n' ...
'8 1 2 0 1 7 8\n' ...
'9 1 2 0 1 8 3\n' ...
'10 1 2 0 2 3 9\n' ...
'11 1 2 0 2 9 10\n' ...
'12 1 2 0 2 10 11\n' ...
'13 1 2 0 2 11 12\n' ...
'14 1 2 0 2 12 13\n' ...
'15 1 2 0 2 13 2\n' ...
'16 2 2 0 4 3 9 16\n' ...
'17 2 2 0 4 2 4 15\n' ...
'18 2 2 0 4 14 17 19\n' ...
'19 2 2 0 4 14 18 20\n' ...
'20 2 2 0 4 5 6 17\n' ...
'21 2 2 0 4 10 11 18\n' ...
'22 2 2 0 4 7 8 19\n' ...
'23 2 2 0 4 12 13 20\n' ...
'24 2 2 0 4 14 16 18\n' ...
'25 2 2 0 4 14 15 17\n' ...
'26 2 2 0 4 14 19 16\n' ...
'27 2 2 0 4 14 20 15\n' ...
'28 2 2 0 4 6 7 25\n' ...
'29 2 2 0 4 11 12 26\n' ...
'30 2 2 0 4 2 21 13\n' ...
'31 2 2 0 4 3 22 8\n' ...
'32 2 2 0 4 9 24 16\n' ...
'33 2 2 0 4 4 23 15\n' ...
'34 2 2 0 4 17 25 19\n' ...
'35 2 2 0 4 18 26 20\n' ...
'36 2 2 0 4 3 16 22\n' ...
'37 2 2 0 4 2 15 21\n' ...
'38 2 2 0 4 9 10 24\n' ...
'39 2 2 0 4 4 5 23\n' ...
'40 2 2 0 4 15 20 21\n' ...
'41 2 2 0 4 16 19 22\n' ...
'42 2 2 0 4 5 17 23\n' ...
'43 2 2 0 4 10 18 24\n' ...
'44 2 2 0 4 8 22 19\n' ...
'45 2 2 0 4 13 21 20\n' ...
'46 2 2 0 4 15 23 17\n' ...
'47 2 2 0 4 16 24 18\n' ...
'48 2 2 0 4 6 25 17\n' ...
'49 2 2 0 4 11 26 18\n' ...
'50 2 2 0 4 12 20 26\n' ...
'51 2 2 0 4 7 19 25\n' ...
'$EndElements\n'];
end

