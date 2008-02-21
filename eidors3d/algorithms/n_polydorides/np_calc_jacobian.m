function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: np_calc_jacobian.m,v 1.14 2008-02-21 19:23:18 aadler Exp $

p= np_fwd_parameters( fwd_model );

s_mat= calc_system_mat( fwd_model, img );

dims= size(fwd_model.nodes,2)
if dims==2
elseif dims==3
   v_f = np_calc_3d_fields( fwd_model, img );

   tol = 1e-5; %tolerance for the forward solver

   % Calculating the Jacobian
   Vfwd = forward_solver(s_mat.E, p.I, tol, s_mat.perm);
   J = jacobian_3d_fields(Vfwd,s_mat.Ela,s_mat.D, p.elec, ...
                          p.vtx,p.simp, img.elem_data, v_f, p.df);
else
   error('dims not 2 or 3');
end

% calculate normalized Jacobian if required
if p.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,p.n_elem));
end

