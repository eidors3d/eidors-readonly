function [elec_nodes, refine_nodes] = dm_mk_elec_nodes( n_elecs, elec_width);
% DM_MK_ELEC_NODES: create node points for dm_mk_fwd_model
% [elec_nodes, refine_nodes] = dm_mk_elec_nodes( n_elecs, elec_width);
%
%  elec_nodes:        cell of matrix N x [x,y,{z}] for each electrode
%  refine_nodes:      vector of fixed nodes to add to model to refine it

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_elec_nodes.m,v 1.1 2008-03-28 19:47:54 aadler Exp $


   radius=1;
   th_delta = elec_width/2/pi/radius;
% electrodes start top and go clockwise
   refine_nodes= [];
   for i=1:n_elecs
      th= (i-1)*2*pi/n_elecs;
      the= th + th_delta*[-1;-.5;0;.5;1];
      this_e_node= radius*[sin(the),cos(the)];
      elec_nodes{i} = this_e_node;

      thr= th + th_delta*[-5;-2;-1;-.5;0;.5;1;2;5];
      csthr= [sin(thr),cos(thr)];
      refine_nodes= [refine_nodes; ...
                   radius*csthr([1:2,8:9],:); ...
               .98*radius*csthr([1:9],:); ...
               .95*radius*csthr([1:2:9],:); ...
               .90*radius*csthr([1:4:9],:)];
   end
