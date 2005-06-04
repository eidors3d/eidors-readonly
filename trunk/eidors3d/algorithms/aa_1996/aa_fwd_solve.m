function data =aa_fwd_solve(fwd_model, img)
% AA_FWD_SOLVE: data= aa_fwd_solve( fwd_model, img)
% Fwd solver for Andy Adler's EIT code
% data = measurements struct
% fwd_model = forward model
% img = image struct

% (C) 1995-2002 Andy Adler
% Ref: Adler & Guardo (1996) IEEE T. Med Imaging
% $Id: aa_fwd_solve.m,v 1.4 2005-06-04 17:49:45 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
s_mat= calc_system_mat( fwd_model, img );

gnd= fwd_model.gnd_node;

d=   pp.n_dims+1;
idx= 1:pp.n_elem*d;
SS= s_mat.SS*sparse(idx,idx, img.elem_data(ceil(idx/d)) );
z= s_mat.CC'* SS * s_mat.CC;

idx= 1:pp.n_node;
idx(gnd) = [];

v= zeros(pp.n_node,pp.n_stim);
v(idx,:)= z(idx,idx) \ pp.QQ(idx,:);
%vv= v(MES(1,:),:)- v(MES(2,:),:);
%vv= vv(ELS);

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
data.name= 'solved by np_fwd_solve';
% TODO: figure out how to describe measurment pattern
data.configuration='unknown';
