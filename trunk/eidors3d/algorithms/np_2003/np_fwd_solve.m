function data= np_fwd_solve( fwd_model, img)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% img = image struct

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: np_fwd_solve.m,v 1.8 2007-04-12 14:54:55 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

s_mat= calc_system_mat( fwd_model, img );

Vfwd = forward_solver(s_mat.E, p.I, tol, s_mat.perm);

Velec=Vfwd( p.n_node+(1:p.n_elec),:);
voltH = zeros( p.n_meas, 1 );
idx=0;
for i=1:p.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   voltH( idx+(1:n_meas) ) = meas_pat*Velec(:,i);
   idx= idx+ n_meas;
end

% create a data structure to return
data.meas= voltH;
data.time= -1; % unknown
data.name= 'solved by np_fwd_solve';
% TODO: figure out how to describe measurment pattern
data.configuration='unknown';
