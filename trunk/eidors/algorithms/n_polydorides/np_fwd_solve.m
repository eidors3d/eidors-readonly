function data= np_fwd_solve( fwd_model, img)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% Input:
%    fwd_model = forward model
%    img       = image struct
% Output:
%    data = measurements struct
% Options:
%    img.fwd_solve.get_all_meas = 1  (return all meas as data.volt)

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

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
try; if img.fwd_solve.get_all_meas == 1
   data.volt = Vfwd(1:p.n_node,:); % but not on CEM nodes
end; end
