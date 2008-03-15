function system_mat = calc_system_mat( fwd_model, img)
% CALC_SYSTEM_MAT: calculate FEM system matrix from fwd_model and image
% 
%    system_mat= calc_system_mat( fwd_model, image)
% OR
%    system_mat= calc_system_mat( image)
%
% it will call the fwd_model.solve
%
% system_mat  
%   system_mat.E    is FEM system_matrix
%   system_mat.perm is permutation of E  i.e. E(perm,perm)
% fwd_model is a fwd_model structure
% image     is an image structure

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_system_mat.m,v 1.17 2008-03-15 22:36:52 aadler Exp $

if nargin>1
   img.fwd_model= fwd_model;
end

system_mat= eidors_obj('get-cache', img, 'system_mat');
if ~isempty(system_mat)
   eidors_msg('system_mat: using cached value', 3);
   return
end

system_mat= feval(img.fwd_model.system_mat, img.fwd_model, img);

eidors_obj('set-cache', img, 'jacobian', system_mat);
eidors_msg('calc_system_mat: setting cached value', 3);

