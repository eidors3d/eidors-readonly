function img= solve_use_matrix( inv_model, data1, data2)
% SOLVE_USE_MATRIX solve using reconstruction matrix
% img= solve_use_matrix( inv_model, data1, data2)
% img        => output image (or vector of images)
% data1      => differential data at earlier time
% data2      => differential data at later time
% inv_model  => inverse model struct
% inv_model.solve_use_matrix.RM  = rec_matrix
% inv_model.solve_use_matrix.map = map to elements (optional)
%    map is the size of the reconst model. contains
%    zero if elem is unfilled, otherwise pointer to matrix row
% inv_model.solve_use_matrix.solve_nodes (default to elements)
%    if 1, then solution goes to the img.node_data value

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

dv = calc_difference_data( data1, data2, inv_model.fwd_model);

RM = inv_model.solve_use_matrix.RM;

if isfield(inv_model.solve_use_matrix,'map')
   map = inv_model.solve_use_matrix.map;
else
   map = 1:size(RM,1);
end


rec= inv_model.solve_use_matrix.RM * dv;
sol = rec(map,:);

% create a data structure to return
try   img.name= ['solved by ', inv_model.name];
catch img.name= ['solved by solve_use_matrix'];
end

solve_to = 'elem_data';
try if inv_model.solve_use_matrix.solve_nodes == 1
   solve_to = 'node_data';
end; end

img.(solve_to) = sol;
img.fwd_model= inv_model.fwd_model;
