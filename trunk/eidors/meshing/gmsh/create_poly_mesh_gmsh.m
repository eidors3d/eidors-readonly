function mdl = create_poly_mesh_gmsh(basename, vertices, elem_size )
% Create a 2D Polygon FEM using GMSH
% mdl= CREATE_POLY_MESH_GMSH(basename, vertices, elem_size)
%
% mdl - EIDORS forward model
% vertices - martix with one row for each vertex
% elem_size - size of the element

% (C) 2011 Bartosz Sawicki. License: GPL version 2 or version 3
% $Id: create_gmsh_2d_circle.m 1535 2008-07-26 15:36:27Z aadler $

geo_filename = sprintf('%s.geo', basename);
geo_fid= fopen(geo_filename,'w');

nv = size(vertices,1);

% Points to define polygon
point_no = 1;
for vi = 1:nv
    fprintf(geo_fid,'Point(%d) = {%f, %f, 0, %f};\n',point_no, ...
        vertices(vi, 1), vertices(vi, 2), elem_size );
    point_no = point_no + 1;
end;

% Edges
line_no = 1;
for vi = 1:(nv-1)
    fprintf(geo_fid,'Line (%d) = {%d,%d};\n', line_no, ...
    vi, vi+1);
    line_no = line_no + 1;
end;
fprintf(geo_fid,'Line(%d) = {%d,%d};\n', line_no, nv,1);
line_no = line_no + 1;

% Polygon as a closed loop
fprintf(geo_fid,'Line Loop(%d) = {', line_no );
for vi = 1:(nv-1)
    fprintf(geo_fid,'%d,', vi);
end;
fprintf(geo_fid,'%d};\n', nv );

line_no = line_no + 1;

fprintf(geo_fid, 'Plane Surface(%d) = {%d};\n', line_no, line_no-1);

fclose(geo_fid);

% Call Gmsh 
disp('Calling Gmsh. Please wait ...');
call_gmsh( geo_filename);

msh_filename = sprintf('%s.msh', basename);

disp(['Now reading back data from file: ' msh_filename])

[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh(msh_filename);
mdl.nodes    = vtx;
mdl.elems    = simp;
mdl.name = 'Polygon mesh';
mdl= eidors_obj('fwd_model', mdl);


