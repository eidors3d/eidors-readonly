function param = np_fwd_parameters( fwd_model )
% NP_FWD_PARAMETERS: data= np_fwd_solve( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Nick Polydorides EIDORS3D code
%   param.n_elem   => number of elements
%   param.n_elec   => number of electrodes
%   param.n_node   => number of nodes (vertices)
%   param.n_stim   => number of current stimulation patterns
%   param.n_meas   => number of measurements (total)
%   param.vtx      => vertex matrix
%   param.simp     => connection matrix
%   param.df       => vector of measurements for each current pattern
%   param.elec     => nodes attached to each electrode
%   param.zc       => vector of contact impedances
%   param.indH     => electrodes used for each measurement
%   param.I        => RHS (current term) for FEM solution
%   param.Ib       => Current for electrodes
%   param.sym      => 'sym' parameter
%   param.gnd_ind  => node attached to ground
% $Id: np_fwd_parameters.m,v 1.1 2004-07-18 03:17:25 aadler Exp $

vtx= fwd_model.nodes;
simp= fwd_model.elems;
% calc num electrodes, nodes, stim_patterns
n_elem= size(simp,1);
n_elec=  length(fwd_model.electrode );
n_node = size(fwd_model.nodes,1);
n_stim = length(fwd_model.stimulation );
n_meas = 0;

% Recreate 'df' from fwd_model.stimulation
df= zeros(n_stim,1);
for i=1:n_stim;
    df(i) = size(fwd_model.stimulation(i).meas_pattern ,1);
    n_meas = n_meas + df(i);
end

elec= zeros(length(fwd_model.electrode ), ...
            length(fwd_model.electrode(1).nodes) );
zc  = zeros(length(fwd_model.electrode ), 1);

for i=1:n_elec
    elec(i,:)= fwd_model.electrode(i).nodes;
    zc(i)    = fwd_model.electrode(i).z_contact;
end

% Recreate 'indH' from fwd_model.stimulation
indH= zeros(n_stim, 2);
idx=0;
for i=1:n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern';
   elpos= rem( find(meas_pat(:)== 1)-1 , n_elec) + 1;
   indH( idx+(1:df(i)),1 ) = elpos;
   elpos= rem( find(meas_pat(:)==-1)-1 , n_elec) + 1;
   indH( idx+(1:df(i)),2 ) = elpos;
   idx= idx+ df(i);
end

% calculate FEM RHS matrix, i.e., the current patterns padded with zeroes 
I = zeros( n_elec + n_node, n_stim );
idx=0;
for i=1:n_stim
   I( n_node + (1:n_elec), i ) = ...
         fwd_model.stimulation(i).stim_pattern;
end
I(fwd_model.gnd_node,:) = 0;
Ib= I( n_node + (1:n_elec), : );

% pack into a parameter return list
param.n_elem   = n_elem;
param.n_elec   = n_elec;
param.n_node   = n_node;
param.n_stim   = n_stim;
param.n_meas   = n_meas;
param.vtx      = vtx;
param.simp     = simp;
param.df       = df;
param.elec     = elec;
param.zc       = zc;
param.indH     = indH;
param.I        = I;
param.Ib       = Ib;
param.sym      = fwd_model.misc.sym;
param.gnd_ind  = fwd_model.gnd_node;
