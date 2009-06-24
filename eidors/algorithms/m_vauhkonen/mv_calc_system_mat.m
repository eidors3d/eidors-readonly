function s_mat= mv_calc_system_mat( fwd_model, img)
% MV_CALC_SYSTEM_MAT: s_mat= mv_calc_system_mat( fwd_model, img)
% System Matrix for Marco Vauhkonen's EIDORS2D code
% s_mat.E   = FEM system matrix
% s_mat.Ela = Normalised volumes of the elements
% s_mat.D   = The sgradients of the shape functions over each element.
% s_mat.Vols= Normalised volums of the elements
% s_mat.perm= permutation of system matrix
% fwd_model = forward model
% img       = image background for system matrix calc

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

s_mat = eidors_obj('get-cache', fwd_model, 'mv_system_mat', img);

if ~isempty(s_mat)
   eidors_msg('mv_calc_system_mat: using cached value', 3);
   return
end

p= mv_fwd_parameters( fwd_model );
sigma = img.elem_data;

[Agrad,Kb,M,S,C]=FemMatrix(p.Node,p.Element,p.z_contact);
if all(all( C ~= p.C ))
   error(['The required measurement pattern is not compatible' ...
         'with that required by eidors2d. Please refer to' ...
         'eidors2d_demo1 for an example']);
end
s_mat= UpdateFemMatrix(Agrad,Kb,M,S,sigma);

eidors_obj('set-cache', fwd_model, 'mv_system_mat', s_mat, img);
eidors_msg('mv_calc_system_mat: setting cached value', 3);
