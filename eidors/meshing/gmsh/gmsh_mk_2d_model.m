function mdl = gmsh_mk_2d_model(varargin)
%GMSH_MK_2D_MODEL create a 2D mesh with GMSH
% mdl = gmsh_mk_2d_model(shape)
%
% SHAPE can be:
%  - xy (Nx2)             : a counter- clockwise list of points in 2D 
%                           defining the outer contour
%  - {xy, xy1, xy2, ...}  : allows specifying additional counter-clockwise 
%                           loops  xy1, xy2, etc, which represent holes in  
%                           the bigger contour xy contour
%  - {..., maxsz}         : specifies the characteristic length at each
%                           point (default: 0.1);
%
% For sufficently large MAXSZ, no new points are created on the boundary.
%
%  See also NG_MK_2D_MODEL

%% (C) 2021 Bartlomiej Grychtol, License: GPL version 2 or version 3
% $Id$

[shape, cl] = process_input(varargin{:});

namestem = tempname;
write_geo_file(shape, cl, namestem);
call_gmsh([namestem '.geo']);
[srf,vtx,fc,bc,simp,edg,mat_ind] = gmsh_read_mesh([namestem '.msh']);
mdl.type     = 'fwd_model';
mdl.name = 'gmsh_2d_model';
mdl.nodes    = vtx;
mdl.elems    = simp;
mdl= eidors_obj('fwd_model', mdl);

delete([namestem '.geo']);
delete([namestem '.msh']);


function fname = write_geo_file(shape, cl, namestem)
fname = [namestem '.geo'];
n_points = 0;
n_lines = 0;
loop_idx = [];
fid= fopen(fname,'w');
for i = 1:length(shape)
    n_sh_pts = size(shape{i},1);
    for p = 1:n_sh_pts
        fprintf(fid,'Point(%d) = {%f, %f, 0, %f};\n',...
            n_points + p, shape{i}(p,1), shape{i}(p,2), cl);
    end
    
    for p = 1:n_sh_pts-1
        fprintf(fid,'Line(%d) = {%d, %d};\n',...
            n_lines + p, n_points + p, n_points + p + 1);
    end
    fprintf(fid,'Line(%d) = {%d, %d};\n',...
        n_lines + n_sh_pts, n_points + n_sh_pts, n_points + 1);
    
    fprintf(fid,'Line Loop (%d) = {%d%s};\n',...
        n_lines + n_sh_pts + 1, n_lines + 1,...
        sprintf(', %d', n_lines + (2:n_sh_pts)));
    
    loop_idx = [loop_idx, n_lines + n_sh_pts + 1];
    n_points = n_points + n_sh_pts;
    n_lines = n_lines + n_sh_pts + 1;
end
str = '';
if numel(loop_idx) > 1
    str = sprintf(', %d', loop_idx(2:end));
end

fprintf(fid,'Plane Surface (%d) = {%d%s};\n',...
        n_lines + 1, loop_idx(1), str);
        
fclose(fid);


function [shape, cl] = process_input(shape)
cl = 0.1;
if ~iscell(shape)
   shape = {shape};
end

if numel(shape) > 1 && numel(shape{end}) == 1
    cl = shape{end};
    shape(end) = [];
end