function data= np_fwd_solve( fwd_model, image)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, image)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% image = image struct
% $Id: np_fwd_solve.m,v 1.4 2004-07-16 18:13:18 aadler Exp $

mat_ref= image.elem_data;

elec= zeros(length(fwd_model.electrode ), ...
            length(fwd_model.electrode(1).nodes) );
zc  = zeros(length(fwd_model.electrode ), 1);

for i=1:length(fwd_model.electrode);
    elec(i,:)= fwd_model.electrode(i).nodes;
    zc(i)    = fwd_model.electrode(i).z_contact;
end

[I,Ib] = set_3d_currents(fwd_model.misc.protocol, ...
                         elec, ...
                         fwd_model.nodes, ...
                         fwd_model.gnd_node, ...
                         fwd_model.misc.no_pl);

%Set the tolerance for the forward solver
tol = 1e-5;

[Eref,D,Ela,ppr] = fem_master_full( ...
                fwd_model.nodes, ...
                fwd_model.elems, ...
                mat_ref, ...
                fwd_model.gnd_node, ...
                elec, ...
                zc, ...
                fwd_model.misc.sym);
[Vfwd] = forward_solver(fwd_model.nodes,Eref,I,tol,ppr);
[voltH,voltV,indH,indV,dfr]=get_3d_meas(elec,fwd_model.nodes,Vfwd,Ib, ...
                fwd_model.misc.no_pl);
dfr = dfr(1:2:length(dfr)); %Taking just the horrizontal measurements

% create a data structure to return
data.meas= voltH;
data.time= 0;
data.name= 'solved by np_fwd_solve';
% TODO: Normally the specification of the measurement sequence
% should be part of the fwd_model, here it is part of the solver
data.misc.indH= indH;
data.misc.df= dfr;
   
