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
% $Id: aa_fwd_parameters.m,v 1.15 2006-08-13 11:15:40 aadler Exp $

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

% calculate element volume and surface area
pp.VOLUME=zeros(e,1);
ones_d = ones(1,d);
ones_d1= ones(1,d-1);
d1fac = prod( 1:d-1 );
d2fac = prod( 1:d-2 );
for i=1:e
    this_elem = pp.NODE(:,pp.ELEM(:,i)); 
    pp.VOLUME(i)= abs(det([ones_d;this_elem])) / d1fac;
end

if isfield(fwd_model,'boundary')
    bdy = fwd_model.boundary;
else
    bdy= dubs3(fwd_model.elems);
end

% Matrix to convert Nodes to Electrodes
N2E = sparse(n_elec, n);
for i=1:n_elec
    elec_nodes = fwd_model.electrode(i).nodes;
%   N2E(i, elec_nodes) = 1/length(elec_nodes);
    srf_area   = get_srf_area( bdy, elec_nodes, fwd_model.nodes);
    N2E(i, elec_nodes) = srf_area/sum(srf_area);
end
  

n_meas= 0; % sum total number of measurements
pp.QQ= sparse(n,p);
for i=1:p
    pp.QQ(:,i) = N2E'* fwd_model.stimulation(i).stim_pattern;
    n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern,1);
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

% get surface area of each element on an electrode
function srf_area   = get_srf_area(bdy, elec_nodes, nodes);
   mbdy= zeros(size(bdy));
   for n= elec_nodes
      mbdy= mbdy + (bdy == n); 
   end 
   e_bdy = find( all(mbdy') );

% get boundary faces which match any node (for pt electrode etc)

   if isempty(e_bdy)
      e_bdy = find( sum(mbdy')>=2 );
   end
   if isempty(e_bdy)
      e_bdy = find( any(mbdy') );
   end

   bdy_els= bdy(e_bdy,:);

   % calculate area of each boundary element
   l_e_bdy= length(e_bdy);
   e_srf_area= zeros(l_e_bdy,1);
   if size(nodes,2) == 3
      for i= 1:l_e_bdy
         pts = nodes(bdy_els(i,:),:);
         e_srf_area(i)= triarea3d(pts);
      end
   elseif size(nodes,2) == 2
      for i= 1:l_e_bdy
         pts = nodes(bdy_els(i,:),:);
         e_srf_area(i)= sqrt(sum(diff(pts).^2));
      end
   else
      error('not 2D or 3d');
   end

   %map boundary element area to each node
   l_elec_nodes = length(elec_nodes);
   srf_area= zeros(1,l_elec_nodes);
   for i= 1:l_elec_nodes
      mbdy= any(bdy_els == elec_nodes(i),2); 
      % 1st order node effect integrates to 1/3 of element area
      srf_area(i) = sum( e_srf_area(mbdy) )/3;
   end
