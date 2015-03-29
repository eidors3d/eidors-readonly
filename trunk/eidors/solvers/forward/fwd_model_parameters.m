function param = fwd_model_parameters( fwd_model )
% FWD_MODEL_PARAMETERS: data= fwd_solve_1st_order( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Andy Adler's EIT code
%   param.n_elem     => number of elements
%   param.n_elec     => number of electrodes
%   param.n_node     => number of nodes (vertices)
%   param.n_stim     => number of current stimulation patterns
%   param.n_elec     => number of electrodes
%   param.n_dims     => dimentions (2= 2D, 3=3D)
%   param.n_meas     => number of measurements (total)
%   param.boundary   => FEM boundary
%   param.NODE       => vertex matrix
%   param.ELEM       => connection matrix
%   param.QQ         => Current into each NODE
%   param.VOLUME     => Volume (or area) of each element
%   param.normalize  => difference measurements normalized?
%   param.N2E        => Node to electrode converter
%
% If the stimulation patterns has a 'interior_sources' field,
%   the node current QQ, is set to this value for this stimulation.

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fwd_model) && strcmp(fwd_model, 'UNIT_TEST'); do_unit_test; return; end

param = eidors_cache(@calc_param,fwd_model,'fwd_model_parameters');


% perform actual parameter calculation
function pp= calc_param( fwd_model );

pp.NODE= fwd_model.nodes';
pp.ELEM= fwd_model.elems';

n= size(pp.NODE,2);        %NODEs
d= size(pp.ELEM,1);        %dimentions+1
e= size(pp.ELEM,2);        %ELEMents
try
   p = length(fwd_model.stimulation );
catch 
   p = 0;
end
try
   n_elec= length( fwd_model.electrode );
catch
   n_elec= 0;
   fwd_model.electrode = [];
end

copt.fstr = 'element_volume';
copt.log_level = 4;
pp.VOLUME= eidors_cache(@element_volume, {pp.NODE, pp.ELEM, e, d}, copt );

if isfield(fwd_model,'boundary')
    bdy = double( fwd_model.boundary ); % double because of stupid matlab bugs
else
    bdy = find_boundary(fwd_model.elems);
end

% Matrix to convert Nodes to Electrodes
% Complete electrode model for all electrodes
%  N2E = sparse(1:n_elec, n+ (1:n_elec), 1, n_elec, n+n_elec);
%  pp.QQ= sparse(n+n_elec,p);
copt.cache_obj = {fwd_model.nodes,fwd_model.elems,fwd_model.electrode};
copt.fstr = 'calculate_N2E';
[N2E,cem_electrodes] = eidors_cache(@calculate_N2E,{fwd_model, bdy, n_elec, n}, copt);

if p>0
  stim = fwd_model.stimulation;
  [pp.QQ, pp.n_meas] = calc_QQ_fast(N2E, stim, p);
end

% pack into a parameter return list
pp.n_elem   = e;
pp.n_elec   = n_elec;
pp.n_node   = n;
pp.n_stim   = p;
pp.n_dims   = d-1;
pp.N2E      = N2E;
pp.boundary = bdy;
pp.normalize = mdl_normalize(fwd_model);


% calculate element volume and surface area
function VOLUME = element_volume( NODE, ELEM, e, d)
   VOLUME=zeros(e,1);
   ones_d = ones(1,d);
   d1fac = prod( 1:d-1 );
   if d > size(NODE,1)
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          VOLUME(i)= abs(det([ones_d;this_elem])) / d1fac;
      end
   elseif d == 3 % 3D nodes in 2D mesh
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          d12= det([ones_d;this_elem([1,2],:)])^2;
          d13= det([ones_d;this_elem([1,3],:)])^2;
          d23= det([ones_d;this_elem([2,3],:)])^2;
          VOLUME(i)= sqrt(d12 + d13 + d23 ) / d1fac;
      end
   elseif d == 2 % 3D nodes in 1D mesh (ie resistor mesh)
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          d12= det([ones_d;this_elem([1],:)])^2;
          d13= det([ones_d;this_elem([2],:)])^2;
          d23= det([ones_d;this_elem([3],:)])^2;
          VOLUME(i)= sqrt(d12 + d13 + d23 ) / d1fac;
      end
   else
      warning('mesh size not understood when calculating volumes')
      VOLUME = NaN;
   end
   % calculate element volume and surface area
   VOLUME=zeros(e,1);
   ones_d = ones(1,d);
   d1fac = prod( 1:d-1 );
   if d > size(NODE,1)
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          VOLUME(i)= abs(det([ones_d;this_elem])) / d1fac;
      end
   elseif d == 3 % 3D nodes in 2D mesh
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          d12= det([ones_d;this_elem([1,2],:)])^2;
          d13= det([ones_d;this_elem([1,3],:)])^2;
          d23= det([ones_d;this_elem([2,3],:)])^2;
          VOLUME(i)= sqrt(d12 + d13 + d23 ) / d1fac;
      end
   elseif d == 2 % 3D nodes in 1D mesh (ie resistor mesh)
      for i=1:e
          this_elem = NODE(:,ELEM(:,i)); 
          d12= det([ones_d;this_elem([1],:)])^2;
          d13= det([ones_d;this_elem([2],:)])^2;
          d23= det([ones_d;this_elem([3],:)])^2;
          VOLUME(i)= sqrt(d12 + d13 + d23 ) / d1fac;
      end
   else
      warning('mesh size not understood when calculating volumes')
      VOLUME = NaN;
   end



function [N2E,cem_electrodes] = calculate_N2E( fwd_model, bdy, n_elec, n);
   cem_electrodes= 0; % num electrodes part of Compl. Elec Model
   N2E = sparse(n_elec, n+n_elec);
   for i=1:n_elec
       try
           elec_nodes = fwd_model.electrode(i).nodes;
       catch
           break; %Not a real electrode so don't include
       end
       if length(elec_nodes) ==1 % point electrode (maybe inside body)
          N2E(i, elec_nodes) = 1;
       elseif length(elec_nodes) ==0
          error('zero length electrode specified');
       else
          bdy_idx= find_electrode_bdy( bdy, [], elec_nodes);

          if ~isempty(bdy_idx) % CEM electrode
             cem_electrodes = cem_electrodes+1;
             N2E(i, n+cem_electrodes) =1;
          else % point electrodes
               % FIXME: make current defs between point electrodes and CEMs compatible
             [bdy_idx,srf_area]= find_electrode_bdy( bdy, ...
                            fwd_model.nodes, elec_nodes);
             N2E(i, elec_nodes) = srf_area/sum(srf_area);
          end
       end
   end
   N2E = N2E(:, 1:(n+cem_electrodes));


function [QQ, n_meas] = calc_QQ_slow(N2E, stim, p)
   QQ = sparse(size(N2E,2),1,p);
   n_meas= 0; % sum total number of measurements
   for i=1:p
       src= zeros(size(N2E,2),1);
       try;  src =        N2E'* stim(i).stim_pattern; end
       try;  src = src +  stim(i).interior_sources;   end
       if all(size(src) == [1,1]) && src==0
          error('no stim_patterns or interior_sources provided for pattern #%d',i);
       end
       
       QQ(:,i) = src;
       n_meas = n_meas + size(stim(i).meas_pattern,1);
   end

function [QQ, n_meas] = calc_QQ_fast(N2E, stim, p)
   QQ = sparse(size(N2E,2),1,p);
   try
   ncols = arrayfun(@(x) size(x.stim_pattern,2), stim);
   end
   if any(ncols>1);
      str = 'multiple columns in stim_pattern for patterns: ';
      error('EIDORS:fwd_model_parameters:stim_pattern', ...
            [str, sprintf('#%d ',find(ncols>1))]);
   end
   idx = 1:p; idx(ncols==0)= [];
   try
   QQ(:,idx) = N2E' * horzcat( stim(:).stim_pattern );
   end

   try
   ncols = arrayfun(@(x) size(x.interior_sources,2), stim);
   end
   if any(ncols>1);
      str = 'multiple columns in interior_sources for patterns: ';
      error('EIDORS:fwd_model_parameters:interior_points',...
            [str, sprintf('#%d ',find(ncols>1))]);
   end
   idx = 1:p; idx(ncols==0)= [];
   try
   QQ(:,idx) = QQ(:,idx) +  N2E' * horzcat( stim(:).interior_sources );
   end

   n_meas = size(vertcat(stim(:).meas_pattern),1);



function do_unit_test
   imdl = mk_common_model('a2c2',16); fmdl = imdl.fwd_model;
   pp = fwd_model_parameters(fmdl);
   [QQ1, n1m] = calc_QQ_slow(pp.N2E, fmdl.stimulation, pp.n_stim);
   [QQ2, n2m] = calc_QQ_fast(pp.N2E, fmdl.stimulation, pp.n_stim);
   unit_test_cmp('calc_QQ', norm(QQ1-QQ2,'fro') + norm(n1m-n2m), 0, 1e-15);

   fmdl.stimulation(8).stim_pattern = fmdl.stimulation(8).stim_pattern*[1,2]; 
   err= 0;
   try;  pp = fwd_model_parameters(fmdl);
   catch e
      if strcmp(e.identifier, 'EIDORS:fwd_model_parameters:stim_pattern');
         err = 1;
      end
   end
   unit_test_cmp('error', err, 1);
