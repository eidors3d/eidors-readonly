function param = aa_fwd_parameters( fwd_model )
% NP_FWD_PARAMETERS: data= np_fwd_solve( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Andy Adler's EIT code
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
% $Id: aa_fwd_parameters.m,v 1.1 2005-06-03 19:09:11 aadler Exp $

param = eidors_obj('get-cache', fwd_model, 'aa_1996_fwd_param');

if ~isempty(param)
   eidors_msg('aa_fwd_parameters: using cached value', 3);
   return
end

param = calc_param( fwd_model );

eidors_obj('set-cache', fwd_model, 'aa_1996_fwd_param', param);
eidors_msg('aa_fwd_parameters: setting cached value', 3);

% perform actual parameter calculation
function pp= calc_param( fwd_model );

pp.NODE= fwd_model.nodes';
pp.ELEM= fwd_model.elems';

d= size(pp.ELEM,1);        %dimentions+1
e= size(pp.ELEM,2);        %ELEMents
p = length(fwd_model.stimulation );

QQ= zeros(e,p);
for i=1:n_stim
    stim= fwd_model.stimulation(i);
    % for each electrode in this stimulation pattern
    for elec= find( stim(i).stim_pattern )
    %   fwd_model.electrode(i).z_contact - unused here
        elec_nodes = fwd_model.electrode(i).nodes;
        l_elec_nodes = length(elec_nodes);

        e_stim = fwd_model.
        -----------------------

        QQ(,i)
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
param.df       = df;
param.elec     = elec;
param.zc       = zc;
param.indH     = indH;
param.I        = I;
param.Ib       = Ib;
param.sym      = fwd_model.misc.sym;
param.gnd_ind  = fwd_model.gnd_node;

% get boundary faces which match nodes
function e_bdy  = bdy_with_nodes(bdy,  elec_nodes );
   mbdy= zeros(size(bdy));
   for n= elec_nodes
      mbdy= mbdy + (bdy == n); 
   end 
   e_bdy = find( all(mbdy') );

