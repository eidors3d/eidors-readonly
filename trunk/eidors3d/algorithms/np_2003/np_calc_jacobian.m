function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: np_calc_jacobian.m,v 1.11 2006-11-15 17:18:23 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

[v_f] = m_3d_fields(p.vtx, p.n_elec, p.indH, ...
                    s_mat.E, tol, p.gnd_ind);

% Calculating the Jacobian
J = jacobian_3d(p.I,p.elec,p.vtx,p.simp,p.gnd_ind, ...
                  img.elem_data, ...
                  p.zc,v_f,p.df,tol, p.perm_sym );

% calculate normalized Jacobian
if p.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,p.n_elem));
end

