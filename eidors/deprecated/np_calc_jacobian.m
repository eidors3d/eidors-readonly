function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

warning('EIDORS:deprecated','NP_CALC_JACOBIAN is deprecated as of 07-Jun-2012. Use CALC_JACOBIAN_ADJOINT instead.');

p= np_fwd_parameters( fwd_model );

s_mat= calc_system_mat( fwd_model, img );

v_f = np_calc_3d_fields( fwd_model, img );

tol = 1e-5; %tolerance for the forward solver

% Calculating the Jacobian
Vfwd = forward_solver(s_mat.E, p.I, tol, s_mat.perm);

if isfield(fwd_model,'coarse2fine');
   J = jacobian_3d_fields(Vfwd,s_mat.Ela,s_mat.D, p.elec, ...
                          p.vtx,p.simp, img.elem_data, v_f, p.df, ...
                          fwd_model.coarse2fine);
   nparam= size(fwd_model.coarse2fine,2);
else 
   J = jacobian_3d_fields(Vfwd,s_mat.Ela,s_mat.D, p.elec, ...
                          p.vtx,p.simp, img.elem_data, v_f, p.df);
   nparam= p.n_elem;
end

% calculate normalized Jacobian if required
if p.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,nparam));
end

