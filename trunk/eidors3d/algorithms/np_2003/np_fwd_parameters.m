function param = np_fwd_parameters( fwd_model )
% NP_FWD_PARAMETERS: data= np_fwd_solve( fwd_model )
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Nick Polydorides EIDORS3D code
%   param.n_elem   => number of elements
%   param.n_elec   => number of electrodes
%   param.n_node   => number of nodes (vertices)
%   param.n_stim   => number of current stimulation patterns
%   param.n_meas   => number of measurements (total)
%   param.vtx      => vertex matrix
%   param.simp     => connection matrix
%   param.srf      => boundary triangles
%   param.df       => vector of measurements for each current pattern
%   param.elec     => nodes attached to each electrode
%   param.zc       => vector of contact impedances
%   param.indH     => electrodes used for each measurement
%   param.I        => RHS (current term) for FEM solution
%   param.Ib       => Current for electrodes
%   param.perm_sym => 'sym' parameter
%   param.gnd_ind  => node attached to ground

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: np_fwd_parameters.m,v 1.11 2005-12-07 13:19:15 aadler Exp $

param = eidors_obj('get-cache', fwd_model, 'np_2003_fwd_param');

if ~isempty(param)
   eidors_msg('np_fwd_parameters: using cached value', 3);
   return
end

param = calc_param( fwd_model );

eidors_obj('set-cache', fwd_model, 'np_2003_fwd_param', param);
eidors_msg('np_fwd_parameters: setting cached value', 3);

% perform actual parameter calculation
function param= calc_param( fwd_model );

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

elec= [];
zc  = zeros(n_elec, 1);

if isfield(fwd_model,'boundary')
    srf = fwd_model.boundary;
else
    srf= dubs3(simp);
end

max_elec_nodes=0;
% get electrode parameters
for i=1:n_elec
    elec_nodes= fwd_model.electrode(i).nodes;
    if length(elec_nodes)>1
       e_bdy  = bdy_with_nodes(srf,  elec_nodes );
       n_bdy  = srf(e_bdy,:)';
    else
       n_bdy= elec_nodes;
    end
% elec is a series of nodes matching bdy faces
    en_list{i}= n_bdy(:)';
    if length(n_bdy) > max_elec_nodes
        max_elec_nodes = length(n_bdy);
    end

% contact impedance
    zc(i)    = fwd_model.electrode(i).z_contact;
end

elec= zeros(n_elec, max_elec_nodes);
for i=1:n_elec
    en= en_list{i};
    elec(i,1:length(en)) = en;
end

% Recreate 'indH' from fwd_model.stimulation
indH= zeros(n_stim, 2);
idx=0;
for i=1:n_stim
   meas_pat= fwd_model.stimulation(i).meas_pattern';

   sourcepos= find(meas_pat(:)== 1);
   sourcepos= rem( sourcepos-1 , n_elec) + 1;

   sinkpos  = find(meas_pat(:)==-1);
   sinkpos  = rem( sinkpos  -1 , n_elec) + 1;

   indH( idx+(1:df(i)) , : ) = [sourcepos, sinkpos];
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
param.srf      = srf;
param.df       = df;
param.elec     = elec;
param.zc       = zc;
param.indH     = indH;
param.I        = I;
param.Ib       = Ib;
param.perm_sym = fwd_model.misc.perm_sym;
param.gnd_ind  = fwd_model.gnd_node;

% get boundary faces which match nodes
function e_bdy  = bdy_with_nodes(bdy,  elec_nodes );
   mbdy= zeros(size(bdy));
   for n= elec_nodes
      mbdy= mbdy + (bdy == n); 
   end 
   e_bdy = find( all(mbdy') );

% get boundary faces which match any node
% Use this for point electrodes where there are no bdy faces
% This is sort of an abuse of the model, but at least it can
% produce a reasonable result for pt electrode mdls.

   if isempty(e_bdy)
      e_bdy = find( sum(mbdy')>=2 );
   end
   if isempty(e_bdy)
      e_bdy = find( any(mbdy') );
   end
