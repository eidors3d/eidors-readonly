function data= np_fwd_solve( fwd_model, img)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, img)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% img = image struct
% $Id: np_fwd_solve.m,v 1.3 2004-07-18 03:17:25 aadler Exp $

p= np_fwd_parameters( fwd_model );

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full( p.vtx, p.simp, ...
                img.elem_data, ...
                p.gnd_ind, p.elec, p.zc, p.sym );

[Vfwd] = forward_solver(p.vtx,Eref,p.I,tol,ppr);

% MODIFIED: the specification of the measurement sequence
% should be part of the fwd_model, here it is part of the solver
%[voltH,voltV,indH,indV,dfr]= ...
%    get_3d_meas(elec,fwd_model.nodes,Vfwd,Ib, ...
%                fwd_model.misc.no_pl);
%
%dfr = dfr(1:2:length(dfr)); %Taking just the horrizontal measurements
%data.misc.indH= indH;
%data.misc.df= dfr;

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
data.time= 0;
data.name= 'solved by np_fwd_solve';
% TODO: figure out how to describe measurment pattern
data.configuration='32 electrodes in 2 planes, adjacent drive';
