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
% $Id: calc_system_mat.m,v 1.15 2007-08-29 09:26:18 aadler Exp $

if strcmp( fwd_model.type , 'image')
    img= fwd_model;
    fwd_model= img.fwd_model;
end

system_mat = eidors_obj('calc-or-cache', fwd_model, ...
                  fwd_model.system_mat, img);

