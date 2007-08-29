function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_calc_jacobian.m,v 1.14 2007-08-29 09:04:05 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

[v_f] = m_3d_fields(p.vtx, p.n_elec, p.indH, ...
                    s_mat.E, tol, p.gnd_ind);

% Calculating the Jacobian
Vfwd = forward_solver(s_mat.E, p.I, tol, s_mat.perm);
J = jacobian_3d_fields(Vfwd,s_mat.Ela,s_mat.D, p.elec, ...
                       p.vtx,p.simp, img.elem_data, v_f, df);

% calculate normalized Jacobian
if p.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,p.n_elem));
end

