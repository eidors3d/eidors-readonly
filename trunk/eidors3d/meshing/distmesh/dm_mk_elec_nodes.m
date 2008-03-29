function [elec_nodes, refine_nodes] = dm_mk_elec_nodes( elec_posn, ...
                    elec_width, refine_level);
% DM_MK_ELEC_NODES: create node points for dm_mk_fwd_model
% [elec_nodes, refine_nodes] = dm_mk_elec_nodes( n_elecs, elec_width);
% INPUT:
%  elec_posn:        vector N x [x,y,{z}] of centre of each electrode
%  elec_width:        vector N x 1 of width of each electrode
%  refine_level:      vector N x 1 of refine_level for each electrode
%      refine_level = 0 -> no refinement
%      refine_level = 1 -> more refinement etc
%
% elec_width and refine_level may be a scalar
%
% OUTPUT:
%  elec_nodes:        cell of matrix N x [x,y,{z}] for each electrode
%  refine_nodes:      vector of fixed nodes to add to model to refine it

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_mk_elec_nodes.m,v 1.2 2008-03-29 00:42:04 aadler Exp $

[ne,nd]= size(elec_posn);
elec_width=   ones(ne,1).*elec_width;   % expand scalar if required
refine_level= ones(ne,1).*refine_level; % expand scalar if required

if nd==2
   [elec_nodes, refine_nodes]= mk_elec_nodes_2d( ...
             elec_posn, elec_width, refine_level, ne); 
elseif nd==3
   error('can`t solve 3D problem, yet');
else
   error('elec_posn isn`t 2 or 3D');
end

function [elec_nodes, refine_nodes]= mk_elec_nodes_2d( ...
             elec_posn, elec_width, refine_level, ne); 

% electrodes start top and go clockwise
   refine_nodes= [];
   for i=1:ne
      [ctr, radius] = find_ctr_rad( elec_posn(i,:), elec_posn);
      th= (i-1)*2*pi/ne;
      th_delta = elec_width(i)/2/pi/radius;
      the= th + th_delta*[-1;-.5;0;.5;1];
      this_e_node= radius*[sin(the),cos(the)];
      this_e_node= this_e_node + ones(size(this_e_node,1),1)*ctr;
      elec_nodes{i} = this_e_node;

      thr= th + th_delta*[-5;-2;-1;-.5;0;.5;1;2;5];
      csthr= [sin(thr),cos(thr)];
      refine_nodes= [refine_nodes; ...
                   radius*csthr([1:2,8:9],:); ...
               .98*radius*csthr([1:9],:); ...
               .95*radius*csthr([1:2:9],:); ...
               .90*radius*csthr([1:4:9],:)];
      refine_nodes= refine_nodes + ones(size(refine_nodes,1),1)*ctr;
   end

function [ctr, rad] = find_ctr_rad( this_elec, elec_posn);
   ctr= [0,0];
   rad= sqrt(sum(this_elec.^2));
