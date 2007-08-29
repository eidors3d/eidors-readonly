function s_mat= np_calc_system_mat( fwd_model, img)
% NP_CALC_SYSTEM_MAT: s_mat= np_calc_system_mat( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% s_mat.E   = FEM system matrix
% s_mat.Ela = Normalised volumes of the elements
% s_mat.D   = The sgradients of the shape functions over each element.
% s_mat.Vols= Normalised volums of the elements
% s_mat.perm= permutation of system matrix
% fwd_model = forward model
% img       = image background for system matrix calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_calc_system_mat.m,v 1.7 2007-08-29 09:04:05 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full( p.vtx, p.simp, ...
                img.elem_data, ...
                p.gnd_ind, p.elec, p.zc, p.perm_sym );

s_mat.E     = Eref;
s_mat.Ela   = Ela;
s_mat.D     = D;
s_mat.Vols  = Ela;
s_mat.perm  = ppr;
