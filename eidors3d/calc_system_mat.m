function system_mat = calc_system_mat( fwd_model, img)
% CALC_SYSTEM_MAT: calculate FEM system matrix from fwd_model and image
% 
%    system_mat= calc_system_mat( fwd_model, image)
%
% it will call the fwd_model.solve
%
% system_mat  
%   system_mat.E    is FEM system_matrix
%   system_mat.perm is permutation of E  i.e. E(perm,perm)
% fwd_model is a fwd_model structure
% image     is an image structure
%
% $Id: calc_system_mat.m,v 1.6 2005-02-23 17:50:32 aadler Exp $

system_mat = eidors_obj('get-cache', fwd_model, 'FEM_system_mat', img);

if ~isempty(system_mat)
   eidors_msg('calc_system_mat: using cached value', 2);
   return
end

system_mat= feval( fwd_model.system_mat, fwd_model, img);

eidors_obj('set-cache', fwd_model, 'FEM_system_mat', system_mat, img);
eidors_msg('calc_system_mat: setting cached value', 2);
