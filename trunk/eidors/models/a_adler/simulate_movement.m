function [vh,vi,xyr_pt]= simulate_movement( img, xyzr );
% SIMULATE_MOVEMENT simulate small conductivity perturbations
% [vh,vi]= simulate_movement( img, xyzr );
%
% Simulate small conductivity perturbations (using the 
%  Jacobian calculation to be accurate and fast).
%
% INPUT:
%   2D Models: specify xyzr = matrix of 3xN.
%      Each perturbation (i) is at (x,y) = xyzr(1:2,i)
%      with radius xyzr(3,i). 
%  
%   3D Models: specify xyzr = matrix of 4xN.
%      Each perturbation (i) is at (x,y,z) = xyzr(1:3,i)
%      with radius xyzr(4,i). 
%
%   xyzr = scalar =N - single spiral of N in medium centre
%  
%   img = eidors image object (with img.fwd_model FEM model).
%
% OUTPUT:
%   vh - homogeneous measurements M x 1
%   vi - target simulations       M x n_points
%   xyr_pt - x y and radius of each point 3 x n_points
%
% Example: Simulate small object in centre:
%    imdl = mk_common_model('b3cr',16);
%    img  = calc_jacobian_bkgnd(imdl);
%    [vh,vi] = simulate_movement( img, [0;0;0;0.05]);

% (C) 2009 Andy Adler. Licensed under GPL v2 or v3
% $Id$

   if size(xyzr) == [1,1]
      path = linspace(0,1,xyzr); phi = 2*pi*path;
      meanodes= mean(    img.fwd_model.nodes  );
      lennodes= size( img.fwd_model.nodes,1); 
      img.fwd_model.nodes = img.fwd_model.nodes - ones(lennodes,1)*meanodes;
      maxnodes= max(max(abs( img.fwd_model.nodes(:,1:2) )));
      img.fwd_model.nodes = img.fwd_model.nodes / maxnodes;

      xyzr = [0.9*path.*sin(phi);0.9*path.*cos(phi);0*path; 0.05 + 0*path];
   end

   Nt = size(xyzr,2);
   c2f = mk_c2f_circ_mapping( img.fwd_model, xyzr);

   img.fwd_model.coarse2fine = c2f;
   J= calc_jacobian( img );

   vh= fwd_solve(img);
   vh=vh.meas;

   vi= vh*ones(1,Nt) + J;

% QUESTON: the accuracy of J will depend on how well we interpolate the
% mesh. How important is this?