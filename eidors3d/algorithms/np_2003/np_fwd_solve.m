function data= np_fwd_solve( fwd_model, image)
% NP_FWD_SOLVE: data= np_fwd_solve( fwd_model, image)
% Fwd solver for Nick Polydorides EIDORS3D code
% data = measurements struct
% fwd_model = forward model
% image = image struct
% $Id: np_fwd_solve.m,v 1.2 2004-07-17 16:41:34 aadler Exp $

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
idx=0;
for i=1:n_stim
   I( n_nodes + (1:n_elec), i ) = ...
         fwd_model.stimulation(i).stim_pattern;
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

% MODIFIED: the specification of the measurement sequence
% should be part of the fwd_model, here it is part of the solver
%[voltH,voltV,indH,indV,dfr]= ...
%    get_3d_meas(elec,fwd_model.nodes,Vfwd,Ib, ...
%                fwd_model.misc.no_pl);
%
%dfr = dfr(1:2:length(dfr)); %Taking just the horrizontal measurements
%data.misc.indH= indH;
%data.misc.df= dfr;

Velec=Vfwd( n_nodes+(1:n_elec),:);
voltH = zeros( n_meas, 1 );
idx=0;
for i=1:n_stim
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
