function J= ms_calc_jacobian( fwd_model, img)
% MS_CALC_JACOBIAN: J= ms_calc_jacobian( fwd_model, img)
% Fwd solver for Manuchehr Soleimani EIDORS3D code
% J         = Jacobian matrix
% fwd_model = forward model
% img = image background for jacobian calc

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: ms_calc_jacobian.m,v 1.3 2007-04-12 14:54:55 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

[v_f] = m_3d_fields(p.vtx, p.n_elec, p.indH, ...
                    s_mat.E, tol, p.gnd_ind);

% calc Vref, or get it from the cache
Vref = eidors_obj('get-cache', fwd_model, 'ms_2005_Vref');
if ~isempty(Vref)
   eidors_msg('ms_calc_jacobian: using cached value', 4);
else
   Vref = forward_solver(s_mat.E,p.I,tol,s_mat.perm);

   eidors_obj('set-cache', fwd_model, 'ms_2005_Vref', Vref);
   eidors_msg('ms_calc_jacobian: setting cached value', 4);
end

% Calculating the Jacobian - use Manuch's calculator
J = jacobian_3d_precalc( ...
                  p.I,p.elec,p.vtx,p.simp,p.gnd_ind, ...
                  Vref, ...
                  p.zc,v_f,p.df,tol, p.perm_sym, ...
                  s_mat.D, s_mat.Vols);

% calculate normalized Jacobian
if p.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,p.n_elem));
end

