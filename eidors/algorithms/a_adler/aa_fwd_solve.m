

function data =aa_fwd_solve(fwd_model, img)
% AA_FWD_SOLVE: data= aa_fwd_solve( fwd_model, img)
% Fwd solver for Andy Adler's EIT code
% data = measurements struct
% fwd_model = forward model
% img = image struct

% (C) 1995-2002 Andy Adler. License: GPL version 2 or version 3
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
% $Id$

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

idx= 1:size(s_mat.E,1);
idx( fwd_model.gnd_node ) = [];

v= zeros(pp.n_node,pp.n_stim);

tol= 1e-5;
v(idx,:)= forward_solver( s_mat.E(idx,idx), pp.QQ(idx,:), tol);

% calc voltage on electrodes
v_els= pp.N2E * v;

% measured voltages from v
vv = zeros( pp.n_meas, 1 );
idx=0;
for i=1:pp.n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   vv( idx+(1:n_meas) ) = meas_pat*v_els(:,i);
   idx= idx+ n_meas;
end


% create a data structure to return
data.meas= vv;
data.time= -1; % unknown
data.name= 'solved by aa_fwd_solve';
data.volt= v(1:pp.n_node,:);
