function data= np_fwd_solve( fwd_model, image)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, image)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% image = image struct
% $Id: np_fwd_solve.m,v 1.6 2004-07-17 14:25:03 aadler Exp $

mat_ref= image.elem_data;

% calc num electrodes, nodes, stim_patterns
n_elec=  length(fwd_model.electrode );
n_nodes= size(fwd_model.nodes,1);
n_stim = length(fwd_model.stimulation );
n_meas = 0;
for i=1:n_stim;
    n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern ,1);
end

elec= zeros(n_elec, length(fwd_model.electrode(1).nodes) );
zc  = zeros(n_elec, 1);

for i=1:n_elec
    elec(i,:)= fwd_model.electrode(i).nodes;
    zc(i)    = fwd_model.electrode(i).z_contact;
end

% calculate FEM RHS matrix, i.e., the current patterns padded with zeroes 
I = zeros( n_elec + n_nodes, n_stim );
M = zeros( n_meas, n_elec );
idx=0;
for i=1:n_stim
   I( n_nodes + (1:n_elec), i ) = ...
         fwd_model.stimulation(i).stim_pattern;
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   M( idx+(1:n_meas),: ) = meas_pat;
   idx= idx+ n_meas;
end
I(fwd_model.gnd_node,:) = 0;
Ib= I( n_nodes + (1:n_elec), : );

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

Velec=Vfwd( n_nodes+(1:n_elec),:);
voltH = zeros( n_meas, 1 );
idx=0;
for i=1:n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern;
   n_meas  = size(meas_pat,1);
   voltH( idx+(1:n_meas) ) = meas_pat*Velec(:,i);
   idx= idx+ n_meas;
end

dfr = dfr(1:2:length(dfr)); %Taking just the horrizontal measurements

% create a data structure to return
data.meas= voltH;
data.time= 0;
data.name= 'solved by np_fwd_solve';
% TODO: Normally the specification of the measurement sequence
% should be part of the fwd_model, here it is part of the solver
data.misc.indH= indH;
data.misc.df= dfr;
   
