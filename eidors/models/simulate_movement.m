function [vh,vi,xyzr,c2f]= simulate_movement( img, xyzr, value );
% SIMULATE_MOVEMENT simulate small conductivity perturbations
% [vh,vi,xyzr, c2f]= simulate_movement( img, xyzr );
%
% Simulate small conductivity perturbations (using the 
%  Jacobian calculation to be fast).
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
%   value = the conductivity perturbation of the target (either scalar, or
%           vector)
%
% OUTPUT:
%   vh - homogeneous measurements M x 1
%   vi - target simulations       M x n_points
%   xyzr - x y and radius of each point 3 x n_points
%   c2f - image representation of the simulations
%
% Example: Simulate small 3D object in centre:
%    img = mk_image( mk_common_model('b3cr',16) ); % 2D Image
%    [vh,vi] = simulate_movement( img, [0;0;0;0.05]);
% Example: Simulate small 2D object in at x near side:
%    img = mk_image( mk_common_model('d2d3c',16) ); % 2D Image
%    [vh,vi] = simulate_movement( img, [0.9;0;0.05]);
%    show_fem(inv_solve(mk_common_model('c2c2',16),vh,vi)) % Reconst (new mesh)
% 

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
   [c2f failed] = mk_c2f_circ_mapping( img.fwd_model, xyzr);
   xyzr(:,failed) = [];
   c2f(:,failed) = [];
   Nt = Nt - nnz(failed);
   img.fwd_model.coarse2fine = c2f;

   % We don't want a normalized jacobian here
   img.fwd_model = mdl_normalize(img.fwd_model,0);

   J= calc_jacobian( img );
   J= move_jacobian_postprocess( J, img, Nt);

   vh= fwd_solve(img);
   vh=vh.meas;

   vi= vh*ones(1,Nt) + J;
   
   % Would this be the slow approach?:
%    vi = vh*zeros(1,Nt);
%    for i = 1: Nt
%        img.elem_data = 1 - c2f(:,i);
%        jnk = fwd_solve(img);
%        vi(:,i) = jnk.meas;
%    end

function J= move_jacobian_postprocess( J, img, Nt)
   if size(J,2) == Nt; % No problem
       return ;
   % Check if movement jacobian introduced electrode movements (elecs * dims)
   elseif size(J,2) == Nt + ...
        length(img.fwd_model.electrode) * size(img.fwd_model.nodes,2)
      J = J(:,1:Nt);
   else
      error('Jacobian calculator is not doing the coarse2fine mapping. This is a bug.');
   end





% QUESTON: the accuracy of J will depend on how well we interpolate the
% mesh. How important is this?
