function data =aa_fwd_solve(fwd_model, img)
% AA_FWD_SOLVE: data= aa_fwd_solve( fwd_model, img)
% Fwd solver for Andy Adler's EIT code
% Input:
%    fwd_model = forward model
%    img       = image struct
% Output:
%    data = measurements struct
% Options: (to return internal FEM information)
%    img.fwd_solve.get_all_meas = 1 (data.volt = all FEM nodes, but not CEM)
%    img.fwd_solve.get_all_nodes= 1 (data.volt = all nodes, including CEM)

% (C) 1995-2002 Andy Adler. License: GPL version 2 or version 3
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
% $Id$

% correct input paralemeters if function was called with only img
if nargin==1 && strcmp(fwd_model.type, 'image');
    img = fwd_model;
    fwd_model= img.fwd_model;
end

pp= fwd_model_parameters( fwd_model );
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
try; if img.fwd_solve.get_all_meas == 1
   data.volt = v(1:pp.n_node,:); % but not on CEM nodes
end; end
try; if img.fwd_solve.get_all_nodes== 1
   data.volt = v;                % all, including CEM nodes
end; end
