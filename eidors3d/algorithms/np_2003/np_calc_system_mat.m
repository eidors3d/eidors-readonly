function s_mat= np_calc_system_mat( fwd_model, img)
% NP_CALC_SYSTEM_MAT: s_mat= np_calc_system_mat( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% s_mat.E   = FEM system matrix
% s_mat.perm= permutation of system matrix
% fwd_model = forward model
% img       = image background for system matrix calc
% $Id: np_calc_system_mat.m,v 1.2 2005-02-23 16:12:18 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full( p.vtx, p.simp, ...
                img.elem_data, ...
                p.gnd_ind, p.elec, p.zc, p.sym );

s_mat.E     = Eref;
s_mat.perm  = ppr;
