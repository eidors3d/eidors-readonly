function param = ls_fwd_model_parameters( fwd_model )
% LS_FWD_MODEL_PARAMETERS: data= fwd_solve_1st_order( fwd_model, image)
% Extract parameters from a 'fwd_model' struct which are 
% appropriate for Andy Adler's EIT code with hex/quad elements
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
% $Id: fwd_model_parameters.m 3931 2013-04-20 15:32:57Z aadler $

param = eidors_obj('get-cache', fwd_model, 'fwd_model_parameters');

if ~isempty(param);
   eidors_msg('fwd_model_parameters: using cached value', 3);
   return
end

param = calc_param( fwd_model );

eidors_obj('set-cache', fwd_model, 'fwd_model_parameters', param);
eidors_msg('fwd_model_parameters: setting cached value', 3);

% perform actual parameter calculation
function pp= calc_param( fwd_model );

pp.NODE= fwd_model.nodes';
pp.ELEM= fwd_model.elems';

n= size(pp.NODE,2);         %NODEs
d0= size(pp.NODE,1);
d1= size(pp.ELEM,1);        %dimentions+1
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
end

% calculate element volume and surface area
pp.VOLUME=zeros(e,1);
for j=1:e;
a=pp.NODE(:,pp.ELEM(:,j));
if d0==2;
        P=[-0.5773 -0.5773;...
            0.5773 -0.5773;...
            0.5773 0.5773;...
            -0.5773 0.5773];
        for i=1:d1
            k=[0 1-P(i,2) P(i,2)-P(i,1) P(i,1)-1;...
               P(i,2)-1 0 P(i,1)+1 -P(i,1)-P(i,2);...
               P(i,1)-P(i,2) -P(i,1)-1 0 P(i,2)+1;...
               1-P(i,1) P(i,2)+P(i,1) -P(i,2)-1 0];
            det_J=(1/8)*a(1,:)*k*a(2,:)';
            pp.VOLUME(j)=pp.VOLUME(j)+det_J;
        end
elseif d0==3;
     P=[-0.5773 -0.5773 0.5773;...
            0.5773 -0.5773 0.5773;...
            0.5773 0.5773 0.5773;...
            -0.5773 0.5773 0.5773;...
            -0.5773 -0.5773 -0.5773;...
            0.5773 -0.5773 -0.5773;...
            0.5773 0.5773 -0.5773;...
            -0.5773 0.5773 -0.5773];
        for i=1:d1;
            dN=(1/8)*[-(1-P(i,2))*(1-P(i,3)) (1-P(i,2))*(1-P(i,3)) (1+P(i,2))*(1-P(i,3)) -(1+P(i,2))*(1-P(i,3))...
                -(1-P(i,2))*(1+P(i,3)) (1-P(i,2))*(1+P(i,3)) (1+P(i,2))*(1+P(i,3)) -(1+P(i,2))*(1+P(i,3));...
                -(1-P(i,1))*(1-P(i,3)) -(1+P(i,1))*(1-P(i,3)) (1+P(i,1))*(1-P(i,3)) (1-P(i,1))*(1-P(i,3))...
                -(1-P(i,1))*(1+P(i,3)) -(1+P(i,1))*(1+P(i,3)) (1+P(i,1))*(1+P(i,3)) (1-P(i,1))*(1+P(i,3));...
                -(1-P(i,1))*(1-P(i,2)) -(1+P(i,1))*(1-P(i,2)) -(1+P(i,1))*(1+P(i,2)) -(1-P(i,1))*(1+P(i,2))...
                (1-P(i,1))*(1-P(i,2)) (1+P(i,1))*(1-P(i,2)) (1+P(i,1))*(1+P(i,2)) (1-P(i,1))*(1+P(i,2))];
            J=dN*a';
            pp.VOLUME(j)=pp.VOLUME(j)+det(J);
        end
else
   warning('mesh size not understood when calculating volumes')
   pp.VOLUME = NaN;
end
end


bdy = ls_find_boundary(fwd_model.elems);


% Matrix to convert Nodes to Electrodes
% Complete electrode model for all electrodes
%  N2E = sparse(1:n_elec, n+ (1:n_elec), 1, n_elec, n+n_elec);
%  pp.QQ= sparse(n+n_elec,p);

cem_electrodes= 0; % num electrodes part of Compl. Elec Model
N2E = sparse(n_elec, n+n_elec);
pp.QQ= sparse(n+n_elec,p);

for i=1:n_elec;
    try
        elec_nodes = fwd_model.electrode(i).nodes;
    catch
        break; %Not a real electrode so don't include
    end
    if length(elec_nodes) ==1; % point electrode (maybe inside body)
       N2E(i, elec_nodes) = 1;
    elseif length(elec_nodes) ==0;
       error('zero length electrode specified');
    else
       bdy_idx= ls_find_electrode_bdy( bdy, [], elec_nodes);

       if ~isempty(bdy_idx); % CEM electrode
          cem_electrodes = cem_electrodes+1;
          N2E(i, n+cem_electrodes) =1;
       else % point electrodes
            % FIXME: make current defs between point electrodes and CEMs compatible
          [bdy_idx,srf_area]= ls_find_electrode_bdy( fwd_model.boundary, ...
                         fwd_model.nodes, elec_nodes);
          N2E(i, elec_nodes) = srf_area/sum(srf_area);
       end
    end
end
N2E = N2E(:, 1:(n+cem_electrodes));
pp.QQ= pp.QQ(1:(n+cem_electrodes),:);


n_meas= 0; % sum total number of measurements

if p>0;
   stim = fwd_model.stimulation;
end

for i=1:p;
    src= 0;
    try;  src = src +  N2E'* stim(i).stim_pattern; end
    try;  src = src +  stim(i).interior_sources;   end
    if all(size(src) == [1,1]) && src==0;
       error('no stim_patterns or interior_sources provided for pattern #%d',i);
    end
    
    pp.QQ(:,i) = src;
    n_meas = n_meas + size(stim(i).meas_pattern,1);
end


% pack into a parameter return list
pp.n_elem   = e;
pp.n_elec   = n_elec;
pp.n_node   = n;
pp.n_stim   = p;
pp.n_dims   = d0;
pp.n_meas   = n_meas;
pp.N2E      = N2E;
pp.boundary = bdy;
pp.normalize = mdl_normalize(fwd_model);
