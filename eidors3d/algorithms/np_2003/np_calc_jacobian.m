function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% img = image background for jacobian calc
% $Id: np_calc_jacobian.m,v 1.3 2004-07-18 03:17:25 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

% HACK: we need a way to cache previous results so that
% things do not need to be recalculated here
[Eref,D,Ela,ppr] = fem_master_full( p.vtx, p.simp, ...
                img.elem_data, ...
                p.gnd_ind, p.elec, p.zc, p.sym );

[v_f] = m_3d_fields(p.vtx,p.n_elec,p.indH,Eref,tol,p.gnd_ind);

% Calculating the Jacobian
J = jacobian_3d(p.I,p.elec,p.vtx,p.simp,p.gnd_ind, ...
                  img.elem_data, ...
                  p.zc,v_f,p.df,tol, p.sym );

