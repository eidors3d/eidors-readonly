function J= np_calc_jacobian( fwd_model, img)
% NP_CALC_JACOBIAN: J= np_calc_jacobian( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc
% $Id: np_calc_jacobian.m,v 1.6 2005-02-23 14:43:46 aadler Exp $

J = eidors_obj('cache', fwd_model, 'jacobian');
if ~isempty(J)
   eidors_msg('np_calc_jacobian: using cached value',1);
   return;
end

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

[v_f] = m_3d_fields(p.vtx, p.n_elec, p.indH, ...
                    s_mat.E, tol, p.gnd_ind);

% Calculating the Jacobian
J = jacobian_3d(p.I,p.elec,p.vtx,p.simp,p.gnd_ind, ...
                  img.elem_data, ...
                  p.zc,v_f,p.df,tol, p.sym );


eidors_obj('cache', fwd_model, 'jacobian', J);
eidors_msg('np_calc_jacobian: setting cached value',1);
