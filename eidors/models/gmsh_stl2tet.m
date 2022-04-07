function mdl = gmsh_stl2tet(stlfile, maxh, extra)
%GMSH_STL2TET creates a tetrahedral mesh from an stl file
% mdl = gmsh_stl2tet(stlfile, maxh, extra) where:
%        mdl - EIDORS fwd_model struct
%    stlfile - either:
%                 - a path to an stl file
%               OR
%                 - a struct with .vertices and .faces 
%       maxh - maximum edge length (default: coarse mesh)
%      extra - extra command line options to gmsh
%
% If stlfile is a struct, stl_write will be called first and an STL file
% written in a temporary location. 
%
% CASHING: Calls are cashed iff stlfile is a struct.
%
% NOTE: Only one surface per file is allowed.
%
% See also CALL_GMSH, STL_WRITE

% (C) Bartlomiej Grychtol, 2012-2021.
% $Id$

if nargin < 3
   extra = [];
end
if nargin > 1 && ~isempty(maxh)
   extra = [' -clmax ' num2str(maxh), ' ', extra];
end

if isstruct(stlfile)
    opt.cache_obj = {stlfile.vertices, stlfile.faces};
    mdl = eidors_cache(@do_gmsh_stl2tet,{stlfile, extra}, opt);
else
    mdl = do_gmsh_stl2tet(stlfile, extra);
end



function mdl = do_gmsh_stl2tet(stlfile, extra)
if isstruct(stlfile)
    stem = tempname;
    stl_write(stlfile, [stem, '.stl'])
    stlname = [stem '.stl'];
else
    stem = strrep(stlfile,'.stl','');
    stlname = stlfile;
end
%TODO: Some of this could be exposed as options (Algorithm, Optimize, ...)
fid = fopen([stem '.geo'],'w');
fprintf(fid,'Merge "%s";\n',stlname);
fprintf(fid,'Surface Loop(1) = {1};\n');
fprintf(fid,'Volume(2) = {1};\n');
fprintf(fid,'Physical Volume(3) = {2};\n');
fprintf(fid,'Mesh.Algorithm3D=4;\n'); %1=delaunay (tetgen) and 4=frontal (netgen)
fprintf(fid,'Mesh.OptimizeNetgen=1;\n');
fclose(fid);

call_gmsh([stem '.geo'], 3, extra);

mdl = gmsh_mk_fwd_model([stem '.msh'],[],[],[]);

delete([stem '.geo']);
delete([stem '.msh']);
if isstruct(stlfile)
    delete(stlname);
end