function system_mat = calc_system_mat( fwd_model, img)
% CALC_SYSTEM_MAT: calculate FEM system matrix from fwd_model and image
% 
% calc_system_mat can be called as
%    system_mat= calc_system_mat( fwd_model, image)
% or
%    system_mat= calc_system_mat( image)
%
% in each case it will call the fwd_model.solve
%                      or image.fwd_model.solve method
%
% system_mat  
%   system_mat.E    is FEM system_matrix
%   system_mat.perm is permutation of E  i.e. E(perm,perm)
% fwd_model is a fwd_model structure
% image     is an image structure
%
% $Id: calc_system_mat.m,v 1.1 2004-07-18 03:52:13 aadler Exp $

if nargin==1
   img= fwd_model;
   fwd_model= img.fwd_model;
end
system_mat= feval( fwd_model.system_mat, fwd_model, img);
