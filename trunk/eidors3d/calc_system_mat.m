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
% $Id: calc_system_mat.m,v 1.7 2005-09-16 03:17:46 aadler Exp $

system_mat = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.system_mat, img);

