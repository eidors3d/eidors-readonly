function system_mat = calc_system_mat( fwd_model, img)
% CALC_SYSTEM_MAT: calculate FEM system matrix from fwd_model and image
% 
%    system_mat= calc_system_mat( fwd_model, image)
% OR
%    system_mat= calc_system_mat( image)
%
% it will call the fwd_model.system_mat
%
% if fwd_model.system_mat is a matrix, calc_system_mat will return this
% matrix
%
% system_mat  
%   system_mat.E    is FEM system_matrix
%   system_mat.perm is permutation of E  i.e. E(perm,perm)
% fwd_model is a fwd_model structure
% image     is an image structure

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$
%system_mat= feval(fwd_model.system_mat, fwd_model, img);return

if nargin == 1
   img= fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling CALC_SYSTEM_MAT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
end
warning off EIDORS:DeprecatedInterface
fwd_model= img.fwd_model;

if isnumeric(fwd_model.system_mat)
   system_mat = fwd_model.system_mat;

else
   
   copt.cache_obj= {fwd_model, img.elem_data};
   copt.fstr = 'system_mat';
   
   try % in case it's a string
       fwd_model.system_mat = str2func(fwd_model.system_mat); 
   end
   
   system_mat = eidors_cache(fwd_model.system_mat,{fwd_model,img},copt);
   
end

warning on EIDORS:DeprecatedInterface
