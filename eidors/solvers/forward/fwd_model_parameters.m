function param = fwd_model_parameters( fwd_model, opt )
% FWD_MODEL_PARAMETERS: data= fwd_solve_1st_order( fwd_model, image)
%   Internal function to extract parameters from a fwd_model
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
%   param.QQ         => Current into each NODE (Neuman Boundary Conditions)
%   param.VV         => Voltage driven into each NODE (Dirichlet BC - only where QQ is NaN)
%   param.YY         => Output Admittance (1/Impedance) of each node current (for each value in QQ)
%   param.VOLUME     => Volume (or area) of each element
%   param.normalize  => difference measurements normalized?
%   param.N2E        => Node to electrode converter
%   param.v2meas     => Convert node voltages to measurements
%
% If the stimulation patterns has a 'interior_sources' field,
%   the node current QQ, is set to this value for this stimulation.

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(fwd_model) && strcmp(fwd_model, 'UNIT_TEST'); do_unit_test; return; end

if nargin < 2
   opt.skip_VOLUME = 0;
else
   assert(ischar(opt),'opt must be a string');
   assert(strcmp(opt,'skip_VOLUME'),'opt can only be ''skip_VOLUME''');
   opt = struct;
   opt.skip_VOLUME = 1;
end

copt.fstr = 'fwd_model_parameters';
copt.log_level = 4;

param = eidors_cache(@calc_param,{fwd_model, opt},copt);


% perform actual parameter calculation
function pp= calc_param( fwd_model, opt )

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

if ~opt.skip_VOLUME
   copt.fstr = 'element_volume';
   copt.log_level = 4;
   pp.VOLUME= eidors_cache(@element_volume, {pp.NODE, pp.ELEM, e, d}, copt );
end

if isfield(fwd_model,'boundary')
    bdy = double( fwd_model.boundary ); % double because of stupid matlab bugs
else
    bdy = find_boundary(fwd_model.elems);
end
try %add system_mat_fields.CEM_boundary if it exists
   bdy = [bdy;fwd_model.system_mat_fields.CEM_boundary];
end

% Matrix to convert Nodes to Electrodes
% Complete electrode model for all electrodes
%  N2E = sparse(1:n_elec, n+ (1:n_elec), 1, n_elec, n+n_elec);
%  pp.QQ= sparse(n+n_elec,p);
copt.cache_obj = {fwd_model.nodes,fwd_model.elems,fwd_model.electrode};
copt.fstr = 'calculate_N2E';
[N2E,cem_electrodes] = eidors_cache(@calculate_N2E,{fwd_model, bdy, n_elec, n}, copt);

stim = fwd_model.stimulation;
if p>0
  [pp.QQ, pp.VV, pp.n_meas] = calc_QQ_fast(N2E, stim, p);
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

pp.v2meas   = get_v2meas(pp.n_elec, pp.n_stim, stim);

function v2meas = get_v2meas(n_elec,n_stim,stim)
    v2meas = sparse(n_elec*n_stim,0);
    for i=1:n_stim
        meas_pat= stim(i).meas_pattern;
        n_meas  = size(meas_pat,1);
        v2meas((i-1)*n_elec + 1: i*n_elec,end+(1:n_meas)) = meas_pat.';
    end
    v2meas = v2meas'; % Conjugate transpose. For derivation
                      % see Chapter 5 of Adler&Holder 2021


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


function [N2E,cem_electrodes] = calculate_N2E( fwd_model, bdy, n_elec, n);
   cem_electrodes= 0; % num electrodes part of Compl. Elec Model
   N2E = sparse(n_elec, n+n_elec);
   for i=1:n_elec
      eleci = fwd_model.electrode(i);
% The faces field is used only if 
      if isfield(eleci,'faces')  && ~isempty(eleci.faces)
        if ~isempty(eleci.nodes)
           eidors_msg('Warning: electrode %d has both faces and nodes',i);
        end
        % This is a CEM electrode
        cem_electrodes = cem_electrodes+1;
        N2Ei= 1;
        N2Ei_nodes = n+cem_electrodes;

      elseif isfield(eleci,'nodes') 
          elec_nodes = fwd_model.electrode(i).nodes;
         [N2Ei,N2Ei_nodes,cem_electrodes] =  ...
             N2Ei_from_nodes(fwd_model, ...
              bdy, elec_nodes, cem_electrodes,n);
      else
          eidors_msg('Warning: electrode %d has no nodes',i);
          break; %Not a real electrode so don't include
      end
      N2E(i, N2Ei_nodes) = N2Ei;
   end
   % Extra nodes are added if nodes for electronics components are added
   try
      extra_nodes = fwd_model.extra_nodes;
   catch
      extra_nodes = 0;
   end
   N2E = N2E(:, 1:(n+cem_electrodes+extra_nodes));

% If N2E can be made a logical (0-1) matrix, do it.
   if all(N2E(find(N2E(:)))==1)
      N2E = logical(N2E);
   end


function [N2Ei,N2Ei_nodes,cem_electrodes]= N2Ei_from_nodes( ...
      fwd_model, bdy, elec_nodes, cem_electrodes,n);
  if length(elec_nodes) ==1 % point electrode (maybe inside body)
     N2Ei = 1;
     N2Ei_nodes = elec_nodes;
  elseif length(elec_nodes) ==0
    error('EIDORS:fwd_model_parameters:electrode','zero length electrode specified');
  else
     bdy_idx= find_electrode_bdy( bdy, [], elec_nodes);

     if ~isempty(bdy_idx) % CEM electrode
        cem_electrodes = cem_electrodes+1;
        N2Ei= 1;
        N2Ei_nodes = n+cem_electrodes;
     else % a set of point electrodes
          % FIXME: make current defs between point electrodes and CEMs compatible
        [bdy_idx,srf_area]= find_electrode_bdy( bdy, ...
                       fwd_model.nodes, elec_nodes);
        s_srf_area =  sum(srf_area);
        if s_srf_area == 0;
           error('Surface area for elec#%d is zero. Is boundary correct?',i);
        end
        N2Ei = srf_area/s_srf_area;
        N2Ei_nodes= elec_nodes;
     end
  end


function [QQ, VV, n_meas] = calc_QQ_slow(N2E, stim, p)
   QQ = sparse(size(N2E,2),p);
   VV = sparse(size(N2E,2),p); N2E0 = N2E>0;
   n_meas= 0; % sum total number of measurements
   for i=1:p
       src= zeros(size(N2E,2),1);
       try;  src =       N2E' * stim(i).stim_pattern; end
       try;  src = src + stim(i).interior_sources;    end
       if all(size(src) == [1,1]) && src==0
          error('no stim_patterns or interior_sources provided for pattern #%d',i);
       end
       
       QQ(:,i) = src;
       n_meas = n_meas + size(stim(i).meas_pattern,1);

       vlt= zeros(size(N2E,2),1);
       try;  vlt =      N2E0' * stim(i).volt_pattern; end
       VV(:,i) = vlt;
   end

function [QQ, VV, n_meas] = calc_QQ_fast(N2E, stim, p)
   try
   ncols = arrayfun(@(x) size(x.stim_pattern,2), stim);
   end
   if any(ncols>1);
      str = 'multiple columns in stim_pattern for patterns: ';
      error('EIDORS:fwd_model_parameters:stim_pattern', ...
            [str, sprintf('#%d ',find(ncols>1))]);
   end
   idx = 1:p; idx(ncols==0)= [];

   QQ = sparse(size(N2E,2),p);
   try
   QQ(:,idx) = N2E' * horzcat( stim(:).stim_pattern );
   end
   VV = sparse(size(N2E,2),p);
   % For voltages, we just need to know which N2E, not the size



   try
   ncols = arrayfun(@(x) size(x.volt_pattern,2), stim);
   end
   if any(ncols>1);
      str = 'multiple columns in volt_pattern for patterns: ';
      error('EIDORS:fwd_model_parameters:volt_pattern', ...
            [str, sprintf('#%d ',find(ncols>1))]);
   end
   idx = 1:p; idx(ncols==0)= [];

   try
   VV(:,idx) = (N2E>0)' * horzcat( stim(:).volt_pattern );
   end

   n_meas = size(vertcat(stim(:).meas_pattern),1);



function do_unit_test
   imdl = mk_common_model('a2c2',16); fmdl = imdl.fwd_model;
   pp = fwd_model_parameters(fmdl);
   [QQ1, VV1, n1m] = calc_QQ_slow(pp.N2E, fmdl.stimulation, pp.n_stim);
   [QQ2, VV2, n2m] = calc_QQ_fast(pp.N2E, fmdl.stimulation, pp.n_stim);
   unit_test_cmp('calc_QQ', norm(QQ1-QQ2,'fro') + norm(n1m-n2m), 0, 1e-15);
   unit_test_cmp('calc_VV1', norm(VV1,'fro'), 0, 1e-15);
   unit_test_cmp('calc_VV2', norm(VV2,'fro'), 0, 1e-15);

   for i=1:6;
      imdl = mk_common_model('a2C0',4); fmdl = imdl.fwd_model;
      switch i
         case 1; fmdl.stimulation(3).stim_pattern = fmdl.stimulation(3).stim_pattern*[1,2]; 
                 expected_err = 'EIDORS:fwd_model_parameters:stim_pattern';
         case 2; fmdl.stimulation(1).stim_pattern = [];
                 expected_err = ''; expected = zeros(45,4);
                 expected(42:45,2:4) = [0,0,1;-1,0,0;1,-1,0;0,1,-1]*10;
                 param = 'QQ';
         case 3; fmdl.electrode(1).nodes = [];
                 expected_err = 'EIDORS:fwd_model_parameters:electrode';
         case 4; fmdl.stimulation(1).volt_pattern = [zeros(3,1);6];
                 expected_err = ''; expected = zeros(45,4); expected(45,1) = 6;
                 param = 'VV';
         case 5; fmdl.stimulation(3).volt_pattern = [ones(4,2)];
                 expected_err = 'EIDORS:fwd_model_parameters:volt_pattern';
         case 6; fmdl.electrode(3).faces = [1,2;2,3;3,1];
                 fmdl.electrode(3).nodes = [];
                 expected_err = ''; expected = zeros(4,45);
                 expected(:,42:45) = eye(4);
                 param = 'N2E';
      end
      err= '';
      try;  pp = fwd_model_parameters(fmdl);
      catch e
         err= e.identifier;
      end
      if length(expected_err)>0;
         unit_test_cmp(['expected error:',num2str(i)], err, expected_err);
      else
         unit_test_cmp(['case:',num2str(i)], full(pp.(param)), expected);
      end
   end
