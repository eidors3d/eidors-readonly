function mdl = gmsh_stl2tet(stlfile, maxh, extra)
%GMSH_STL2TET creates a tetrahedral mesh from an stl file
% mdl = gmsh_stl2tet(stlfile, maxh, extra) where:
%        mdl - EIDORS fwd_model struct
%    stlfile - path to stl file
%       maxh - maximum edge length (default: coarse mesh)
%      extra - extra command line options to gmsh
%
% NOTE: Only one surface per file is allowed.
%
% See also CALL_GMSH

% (C) Bartlomiej Grychtol, 2012.
% $Id$

if nargin < 3
   extra = [];
end
if nargin > 1 && ~isempty(maxh)
   extra = [' -clmax ' num2str(maxh), ' ', extra];
end

stem = strrep(stlfile,'.stl','');
%TODO: Some of this could be exposed as options (Algorithm, Optimize, ...)
fid = fopen([stem '.geo'],'w');
fprintf(fid,'Merge "%s";\n',stlfile);
fprintf(fid,'Surface Loop(1) = {1};\n');
fprintf(fid,'Volume(1) = {1};\n');
fprintf(fid,'Physical Volume(''object'') = {1};\n');
fprintf(fid,'Mesh.Algorithm3D=4;\n'); %1=delaunay (tetgen) and 4=frontal (netgen)
fprintf(fid,'Mesh.OptimizeNetgen=1;\n');
fclose(fid);

call_gmsh([stem '.geo'], 3, extra);

mdl = gmsh_mk_fwd_model([stem '.msh'],[],[],[]);