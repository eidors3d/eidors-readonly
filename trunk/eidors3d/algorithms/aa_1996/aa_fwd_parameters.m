function param = aa_fwd_parameters( fwd_model )
% AA_FWD_PARAMETERS: data= aa_fwd_solve( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Andy Adler's EIT code
%   param.n_elem     => number of elements
%   param.n_elec     => number of electrodes
%   param.n_node     => number of nodes (vertices)
%   param.n_stim     => number of current stimulation patterns
%   param.n_dims     => dimentions (2= 2D, 3=3D)
%   param.n_meas     => number of measurements (total)
%   param.NODE       => vertex matrix
%   param.ELEM       => connection matrix
%   param.QQ         => Current into each NODE
%   param.VOLUME     => Volume (or area) of each element
%   param.normalize  => difference measurements normalized?

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: aa_fwd_parameters.m,v 1.13 2005-10-27 13:28:08 aadler Exp $

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

n= size(pp.NODE,2);        %NODEs
d= size(pp.ELEM,1);        %dimentions+1
e= size(pp.ELEM,2);        %ELEMents
p = length(fwd_model.stimulation );
n_elec= length( fwd_model.electrode );

% Matrix to convert Nodes to Electrodes
N2E = sparse(n_elec, n);
for i=1:n_elec
    elec_nodes = fwd_model.electrode(i).nodes;
    N2E(i, elec_nodes) = 1/length(elec_nodes);
end
  

n_meas= 0; % sum total number of measurements
pp.QQ= sparse(n,p);
for i=1:p
    pp.QQ(:,i) = N2E'* fwd_model.stimulation(i).stim_pattern;
    n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern,1);
end

% calculate element volume
pp.VOLUME=zeros(e,1);
ones_d = ones(1,d);
d1fac = prod( 1:d-1 );
for i=1:e
    this_elem = pp.NODE(:,pp.ELEM(:,i)); 
    pp.VOLUME(i)= abs(det([ones_d;this_elem])) / d1fac;
end


% pack into a parameter return list
pp.n_elem   = e;
pp.n_elec   = n_elec;
pp.n_node   = n;
pp.n_stim   = p;
pp.n_dims   = d-1;
pp.n_meas   = n_meas;
pp.N2E      = N2E;

if isfield(fwd_model,'normalize_measurements')
   pp.normalize = fwd_model.normalize_measurements;
else
   pp.normalize = 0;
end
