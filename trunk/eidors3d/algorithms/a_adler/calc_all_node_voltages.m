function node_v= calc_all_nodes_voltages( himg );
% CALC_ALL_NODES_VOLTAGES - calculate voltage on all nodes
% node_v= calc_all_nodes_voltages( himg );
% img        => image object
%
% node_v     = n_nodes x n_stims voltage on each node

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_all_node_voltages.m,v 1.2 2008-03-10 18:27:25 aadler Exp $

   % create one "measurement electrode" per node and
   % stim patterns to match
   n_nodes = size(himg.fwd_model.nodes,1);
   n_elecs = length(himg.fwd_model.electrode);
   ze= zeros(n_elecs,1);
   for i= 1:n_nodes
      himg.fwd_model.electrode(n_elecs + i).nodes= i;
      himg.fwd_model.electrode(n_elecs + i).z_contact= 0;
   end

   n_stims = length(himg.fwd_model.stimulation);
   zn= zeros(n_nodes,1);
   zs= zeros(n_stims,1);
   meas_pat = spdiags(zn+1, n_elecs, n_nodes, n_nodes+n_stims);
   for i= 1:n_stims
      himg.fwd_model.stimulation(i).meas_pattern= meas_pat;
      himg.fwd_model.stimulation(i).stim_pattern= ...
         [himg.fwd_model.stimulation(i).stim_pattern; zn];
   end

   node_v= fwd_solve( himg );
   node_v= reshape(node_v.meas,n_nodes,n_stims);



