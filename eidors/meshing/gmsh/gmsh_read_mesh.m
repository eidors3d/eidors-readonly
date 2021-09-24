function [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names,entities] = gmsh_read_mesh(filename)
%[srf,vtx,fc,bc,simp,edg,mat_ind,phys_names,entities] = gmsh_read_mesh(filename)
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
%
% This mostly works on GMSH v2. A very basic GMSH v4 reader is now
% included.

% $Id$
% (C) 2009 Bartosz Sawicki. Licensed under GPL V2
% Modified by James Snyder, Mark Campbell, Symon Stowe, Alistair Boyle

if ischar(filename) && strcmp(filename,'UNIT_TEST'); do_unit_test; return; end

fid = fopen(filename,'r');
phys_names = [];
entities = [];
while 1
    tline = fgetl(fid);
    if ~ischar(tline); fclose(fid); break; end

    if strcmp(tline,'$MeshFormat')
       gmshformat = parse_format( fid );
    elseif strcmp(tline, '$Entities')
       entities= parse_entities( fid, gmshformat );
    elseif strcmp(tline,'$Elements')
       elements= parse_elements( fid, gmshformat );
    elseif strcmp(tline,'$Nodes')
       nodes= get_lines_with_nodes( fid, gmshformat );
    elseif strcmp(tline,'$PhysicalNames')
       phys_names= parse_names( fid, gmshformat );
    end
end

% look up nodes for each of phys_names
if ~isempty(phys_names)
    if (gmshformat >= 4.0)
        assert(~isempty(entities), 'expected $Entities section (GMSH format 4)');
        for i = 1:length(phys_names)
            phys_tag = phys_names(i).tag;
            dim = phys_names(i).dim;
            tmpentities = find(arrayfun(@(x) any(x.phys_tag == phys_tag), entities));
            tags = cat(1,entities(tmpentities).entity_tag);
            tmpelements = find(arrayfun(@(x) and(any(x.entity_tag == tags), (x.dim == dim)), elements));
            phys_names(i).nodes = cat(1,elements(tmpelements).simp);
            for j = tmpelements(:)'
                assert(length(elements(j).phys_tag) > 0, ...
                    'GMSH format v4 mesh volume elements can only have one $PhysicalName when imported into EIDORS');
                elements(j).phys_tag = phys_tag;
            end
        end
    else % GMSH 2.0 doesn't have Entities
        assert(isempty(entities), 'unexpected $Entities section (GMSH format 2)');
        for i = 1:length(phys_names)
            dim = phys_names(i).dim;
            tag = phys_names(i).tag;
            tmpelements = find(arrayfun(@(x) and(any(x.phys_tag == tag), (x.dim == dim)), elements));
            phys_names(i).nodes = cat(1,elements(tmpelements).simp);
        end
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
    mat_ind = cat(1,elements(tet).phys_tag);
else
    vtx = nodes(:,2:3);
    tri = find(arrayfun(@(x)x.type==2,elements));
    simp = cat(1,elements(tri).simp);
    srf = [];
    mat_ind = cat(1,elements(tri).phys_tag);
end

elemtags = cat(1,elements.phys_tag);
fc = elemtags(tri,1);
end

function mat = get_lines_with_nodes( fid, gmshformat )
	tline = fgetl(fid);
	n_rows = parse_rows(tline,gmshformat);
    switch floor(gmshformat)
    % Version 2 Line Format:
    % node-number x-coord y-coord z-coord
    % Version 4 Line Format: (not always like this)
    % node-number
    % x-coord y-coord z-coord
	case 2; mat= fscanf(fid,'%f',[4,n_rows])';
	case 4;
        n = sscanf(tline, '%d')';
        n_block = n(1);
        n_nodes = n(2);
        mat = zeros(n_nodes,4);
        while n_block > 0
            n_block = n_block - 1;
            tline = fgetl(fid);
            blk = sscanf(tline, '%d')';
            blk_nodes = blk(end); % n(2) for v4.0, n(4) for v4.1
            node_tag = zeros(1, blk_nodes);
            el = zeros(blk_nodes,3);
            if gmshformat == 4.1 % v4.1: node tags first as a block, then node coordinates in a block
                for i = 1:blk_nodes
                    tline = fgetl(fid);
                    node_tag(i) = sscanf(tline, '%d')';
                end
                for i = 1:blk_nodes
                    tline = fgetl(fid);
                    el(i,:) = sscanf(tline, '%f')';
                end
            else % v4.0: node tag and coordinates on the same line
                for i = 1:blk_nodes
                    tline = fgetl(fid);
                    data = sscanf(tline, '%f')';
                    node_tag(i) = data(1);
                    el(i,:) = data(2:end);
                end
                %data = fscanf(fid,'%f',[4,blk_nodes])';
                %node_tag = data(:,1);
                %el = data(:,2:end);
            end
            mat(node_tag,:) = [node_tag(:), el(:,1:end)];
        end
        assert(n_block == 0, 'failed to consume all $Nodes blocks');
    otherwise; error('cant parse gmsh file of this format');
	end
end

function gmshformat = parse_format(fid)
   tline = fgetl(fid);
   rawformat = sscanf(tline,'%f');
   tline = fgetl(fid); % should be EndMeshFormat
   gmshformat = rawformat(1);
end

function n_rows = parse_rows(tline, gmshformat)
   n_rows = sscanf(tline,'%d');
   switch floor(gmshformat)
     case 2; n_rows = n_rows(1);
     case 4; n_rows = n_rows(2);
     otherwise; error('cant parse gmsh file of this format');
   end
end

function names = parse_names( fid, version )
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
end % end function

function elements = parse_entities( fid, gmshformat )
    elements = [];
    switch floor(gmshformat)
      case 2; warning('ignoring $Entities$ for GMSH v2 file');
      case 4;
          elements = parse_v4_entities(fid, gmshformat);
      otherwise error('cant parse this file type');
    end
 end

function entities = parse_v4_entities(fid, gmshformat)
% http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
% $Entities
%  numPoints(size_t) numCurves(size_t)
%    numSurfaces(size_t) numVolumes(size_t)
%  pointTag(int) X(double) Y(double) Z(double)
%    numPhysicalTags(size_t) physicalTag(int) ...
%  ...
%  curveTag(int) minX(double) minY(double) minZ(double)
%    maxX(double) maxY(double) maxZ(double)
%    numPhysicalTags(size_t) physicalTag(int) ...
%    numBoundingPoints(size_t) pointTag(int) ...
%  ...
%  surfaceTag(int) minX(double) minY(double) minZ(double)
%    maxX(double) maxY(double) maxZ(double)
%    numPhysicalTags(size_t) physicalTag(int) ...
%    numBoundingCurves(size_t) curveTag(int) ...
%  ...
%  volumeTag(int) minX(double) minY(double) minZ(double)
%    maxX(double) maxY(double) maxZ(double)
%    numPhysicalTags(size_t) physicalTag(int) ...
%    numBoundngSurfaces(size_t) surfaceTag(int) ...
%  ...
% $EndEntities
% Entities are necessary to map '$PhysicalNames' to 'entityDim' & 'entityTag'
% in $Elements and $Nodes.
% Note that the entityTag can repeat for each type of entityDim, so when we
% look up the elements associated with a PhysicalName we need to then use the
% entityDim AND entityTag to find matching Nodes/Elements.
    entities = struct('phys_tag',{},'dim',{},'entity_tag',{});
    tline = fgetl(fid);
    bl = sscanf(tline, '%d')'; % Get the line info
    n_points = bl(1);
    n_curves = bl(2);
    n_surf = bl(3);
    n_vol = bl(4);
    %fprintf('entities: %d pts, %d curves, %d surf, %d vol\n', n_points, n_curves, n_surf, n_vol);
    tline = fgetl(fid); % get the EndElements
    while ~strcmp(tline, '$EndEntities')
        bl = sscanf(tline, '%f')'; % Get the line info
        off = 8; % offset for 'numPhysicalTags'
        if n_points > 0
            n_points = n_points - 1;
            dim = 0; % point
            off = 5; % offset for 'numPhysicalTags'
        elseif n_curves > 0
            n_curves = n_curves - 1;
            dim = 1; % line/curve
        elseif n_surf > 0
            n_surf = n_surf - 1;
            dim = 2; % surface
        elseif n_vol > 0
            n_vol = n_vol - 1;
            dim = 3; % volume
        else
            error('extra $Entities found in v4.1 GMSH file');
        end
        if bl(off) > 0
            entities(end+1).phys_tag = bl((off+1):int64(off+bl(off)));
            entities(end).entity_tag = int64(bl(1));
            entities(end).dim = dim;
        end
        tline = fgetl(fid);
    end
    assert(n_points == 0, 'missing Entity/points from GMSH 4.1 file');
    assert(n_curves == 0, 'missing Entity/curves from GMSH 4.1 file');
    assert(n_surf == 0, 'missing Entity/surfaces from GMSH 4.1 file');
    assert(n_vol == 0, 'missing Entity/volumes from GMSH 4.1 file');
end


function elements = parse_elements( fid, gmshformat )
   tline = fgetl(fid);
   n_rows = parse_rows(tline,gmshformat);
   switch floor(gmshformat)
     case 2; elements = parse_v2_elements(fid,n_rows);
     case 4; elements = parse_v4_elements(fid,tline,gmshformat);
     otherwise error('cant parse this file type');
   end
end

function elements = parse_v4_elements(fid,tline,gmshformat)
% http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format
% $Elements
% numEntityBlocks numElements minElementTag maxElementTag
% ...
% entityDim entityTag elementType numElementsInBlock
% elementTag nodeTag ... nodeTag
% ...
% $EndElements
% An entity is a node, curve, surface or volume
% Each entity has a list of contained elements and corresponding node tags
% 0D - 1 node tag (ignore these)
% 1D - 2 node tags
% 2D - 3 node tags
% 3D - 4 node tags
    elements = struct('simp',{},'type',{},'entity_tag',{},'dim',{},'phys_tag',{});
    bl = sscanf(tline, '%d')'; % Get the line info
    n_blocks = bl(1);
    n_elems = bl(2);
    e_block = 0;
    tline = fgetl(fid);
    while ~strcmp(tline, '$EndElements')
        bl = sscanf(tline, '%d')'; % Get the line info
        if e_block == 0
            n_blocks = n_blocks - 1;
            if gmshformat >= 4.1
                e_dim = bl(1); % entityDim
                e_tag = bl(2); % entityTag
            else % 4.0
                e_tag = bl(1); % entityTag
                e_dim = bl(2); % entityDim
            end
            e_type = bl(3); % elementType
            e_block = bl(4); % numElementsInBlock: Size of entity block
        else
            n_elems = n_elems - 1;
            e_block = e_block - 1;
            elements(bl(1)).type = e_type;
            elements(bl(1)).entity_tag = e_tag;
            elements(bl(1)).dim = e_dim;
            elements(bl(1)).phys_tag = 0; % determined later
            elements(bl(1)).simp = bl(2:end);
        end
        tline = fgetl(fid);
    end
    tline = fgetl(fid); % get the EndElements
    assert(n_elems == 0, 'missing Elements from GMSH 4 file');
    assert(n_blocks == 0, 'missing Element/blocks from GMSH 4 file');
end

function elements = parse_v2_elements(fid,n_rows)
% Line Format:
% elm-number elm-type number-of-tags < tag > ... node-number-list
elements(n_rows).simp = [];
elements(n_rows).phys_tag = [];
elements(n_rows).geom_tag = [];
elements(n_rows).type = [];
elements(n_rows).dim = [];

for i = 1:n_rows
    tline = fgetl(fid);
    n = sscanf(tline, '%d')';
    elements(i).simp = n(n(3) + 4:end);
    elements(i).type = n(2);
    if n(3) > 0 % get tags if they exist
        % By default, first tag is number of parent physical entity
        % second is parent elementary geometrical entity
        % third is number of parent mesh partitions followed by
        % partition ids
        tags = n(4:3+n(3));
        if length(tags) >= 1
            elements(i).phys_tag = tags(1);
            elements(i).dim = length(elements(i).simp) - 1;
            if length(tags) >= 2
                elements(i).geom_tag = tags(2);
            end
        end
    end
end
end

function do_unit_test
   selfdir = fileparts(which('gmsh_read_mesh'));

   [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(fullfile(selfdir, 'test-4.0.msh'));
    unit_test_cmp('v4 vtx ',vtx(2:3,:),[1,0;-1,0])
    unit_test_cmp('v4 simp',simp(2:3,:),[3,7,12; 12, 7,14]);

   [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = gmsh_read_mesh(fullfile(selfdir, 'test-2.2.msh'));
    unit_test_cmp('v2 vtx ',vtx(2:3,:),[1,0;-1,0])
    unit_test_cmp('v2 simp',simp(2:3,:),[2,4,15; 14,17,19]);

   vers = {'2.2', '4.0', '4.1'};
   for ver = vers(:)'
       ver = ver{1};
       [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = ...
           gmsh_read_mesh( fullfile(selfdir, ['box-' ver '.msh']) );
       unit_test_cmp(['2d v' ver ' vtx '],vtx,[0,0;1,0;1,1;0,1;0.5,0.5])
       unit_test_cmp(['2d v' ver ' simp'],simp,[2,5,1;1,5,4;3,5,2;,4,5,3]);
       unit_test_cmp(['2d v' ver ' phys'],{phys_names(:).name},{'elec#1','elec#2','main'});

       [srf,vtx,fc,bc,simp,edg,mat_ind,phys_names] = ...
           gmsh_read_mesh( fullfile(selfdir, ['cube-' ver '.msh']) );
       unit_test_cmp(['3d v' ver ' vtx '],vtx([1,2,13,end],:),[0,0,1;0,0,0;,0.5,0.5,0;0.5,0.5,1])
       unit_test_cmp(['3d v' ver ' simp'],simp([1,2,23,end],:), ...
           [10,11,12,13;9,12,14,11;12,14,10,7;13,10,11,6]);
       unit_test_cmp(['3d v' ver ' phys'],{phys_names(:).name},{'elec#1','elec#2','main'});
   end

end
