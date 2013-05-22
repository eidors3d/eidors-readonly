function mdl = stl_read(fname)
%STL_READ  Read in an stl file and output an EIDORS model struct
% This function is a wrapper for the Matlab Central contribution by 
% Adam H. Aitkenhead called READ_stl from "Mesh voxelization"
% http://www.mathworks.de/matlabcentral/fileexchange/27390-mesh-voxelisation
%
% See also stl_write

% Copyright 2013 Bartłomiej Grychtol
% License: GPL version 2 or 3
% $Id$

if isstr(fname) && strcmp(fname, 'UNIT_TEST'), do_unit_test; return, end

[V] = READ_stl(fname);
N = reshape(shiftdim(V,1),3,[])';
[nodes m n] = unique(N,'rows');
mdl.nodes = nodes;
mdl.elems = reshape(n',3,[])';
mdl.type  = 'fwd_model'; % sort of
mdl.boundary = mdl.elems;


function do_unit_test
   a = [
   -0.8981   -0.7492   -0.2146    0.3162    0.7935    0.9615    0.6751    0.0565   -0.3635   -0.9745
    0.1404    0.5146    0.3504    0.5069    0.2702   -0.2339   -0.8677   -0.6997   -0.8563   -0.4668 ]';
fmdl = ng_mk_extruded_model({2,{a,0.5*a,0.2*a},1},[],[]);

fv.faces    = find_boundary(fmdl);
fv.vertices = fmdl.nodes;
stl_write(fv, 'test.stl');
mdl = stl_read('test.stl');
show_fem(mdl)
delete('test.stl');


