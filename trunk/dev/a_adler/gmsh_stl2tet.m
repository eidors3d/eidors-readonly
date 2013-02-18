function mdl = gmsh_stl2tet(stlfile)
%GMSH_STL2TET creates a tetrahedral mesh from an stl file
%  Only one surface per file is allowed.

% (C) Bartlomiej Grychtol, 2012.
% $Id$

stem = strrep(stlfile,'.stl','');

fid = fopen([stem '.geo'],'w');
fprintf(fid,'Merge "%s";\n',stlfile);
fprintf(fid,'Surface Loop(1) = {1};\n');
fprintf(fid,'Volume(1) = {1};\n');
fprintf(fid,'Physical Volume(''object'') = {1};\n');
fprintf(fid,'Mesh.Algorithm3D=4;\n'); %1=delaunay (tetgen) and 4=frontal (netgen)
fprintf(fid,'Mesh.Optimize=1;\n');
fclose(fid);

call_gmsh([stem '.geo'], 3);

mdl = gmsh_mk_fwd_model([stem '.msh'],[],[],[]);