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
% $Id$

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
      [ctr, radius] = find_ctr_rad( i, elec_posn);
      th= (i-1)*2*pi/ne;
      th_delta = elec_width(i)/2/pi/radius;
      switch refine_level(i)
        case 0,
           the= th;
           thr= th;
        case 1,
           the= th + th_delta*[-1;;0;1];
           thr= th + th_delta*[-3;-1;0;1;3];
        case 2,
           the= th + th_delta*[-1;;0;1];
           thr= th + th_delta*[-5;-2;-1;0;1;2;5];
        case 3,
           the= th + th_delta*[-1;-.5;0;.5;1];
           thr= th + th_delta*[-5;-2;-1;-.5;0;.5;1;2;5];
        case 4,
           the= th + th_delta*[-1;-.5;0;.5;1];
           thr= th + th_delta*[-5;-3;-2;-1.5;-1;-.5;0;.5;1;1.5;2;3;5];
        otherwise
           error('refine level = %d not understood', refine_level(i));
      end
      this_e_node= radius*[sin(the),cos(the)];
      this_e_node= this_e_node + ones(size(this_e_node,1),1)*ctr;
      elec_nodes{i} = this_e_node;

      csthr= [sin(thr),cos(thr)];
      switch refine_level(i)      
        case 0,
           % no refine_nodes
           refine_new= zeros(0,2);
        case 1,
           refine_new= [radius*csthr([1,5],:); ...
                    .97*radius*csthr([3],:)];
        case 2,
           refine_new= [radius*csthr([1:2,6:7],:); ...
                    .98*radius*csthr([1:7],:); ...
                    .95*radius*csthr([1:3:7],:)];
        case 3,
           refine_new= [radius*csthr([1:2,8:9],:); ...
                    .98*radius*csthr([1:9],:); ...
                    .95*radius*csthr([1:2:9],:); ...
                    .90*radius*csthr([1:4:9],:)];
        case 4,
           refine_new= [radius*csthr([1:4,10:13],:); ...
                    .98*radius*csthr([1:13],:); ...
                    .96*radius*csthr([1:2:13],:); ...
                    .93*radius*csthr([1,4,6,8,10,13],:); ...
                    .90*radius*csthr([2:5:12],:)];
        otherwise
           error('refine level = %d not understood', refine_level(i));
      end
      refine_new = refine_new + ones(size(refine_new,1),1)*ctr;
      refine_nodes= [refine_nodes; refine_new];
   end

% Find the ctr and radius of this_elec and the two closest
%  ones. This will allow fitting the electrode to the curvature
%  locally.
% Alg: http://www.geocities.com/kiranisingh/center.html
function [ctr, rad] = find_ctr_rad( idx, elec_posn);
   nelec= size(elec_posn,1);
   idx = rem(idx+[-1,0,1]+nelec-1,nelec)+1;
   x = elec_posn(idx,1);
   y = elec_posn(idx,2);
   s21= x(2)^2 + y(2)^2 - x(1)^2 - y(1)^2;
   s31= x(3)^2 + y(3)^2 - x(1)^2 - y(1)^2;
   x21 = x(2) - x(1);
   x31 = x(3) - x(1);
   y21 = y(2) - y(1);
   y31 = y(3) - y(1);
   n1 = det([s21,y21;s31,y31]);
   n2 = det([x21,s21;x31,s31]);
   D  = det([x21,y21;x31,y31]);
   ctr= [n1,n2]/2/D;
   rad= sqrt((x-ctr(1)).^2 + (y-ctr(2)).^2);
   if std(rad)/mean(rad)>.001; error('PROBLEM WITH ALGORITHM');end
   rad=mean(rad);
