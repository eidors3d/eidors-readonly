function mdl_pts = interp_mesh( mdl, n_interp)
% INTERP_MESH: calculate interpolation points onto mdl elements
%    mdl_pts = interp_mesh( fwd_model, n_interp)
% INPUT:
%    fwd_model: fwd_model structure
%    n_interp:  order of interpolation
%      n_interp = 0 - output elem centres (default)
%      n_interp >=1 - output multiple points per elem
%           in 2D: (n_int+1)*(n_int+2)/2 points per elem
%           in 3D: (n_int+1)*(n_int+2)*(n_int+3)/6 points per elem
% OUTPUT:
%    mdl_pts = N_elems x N_dims x N_points
%   example: for mdl_pts= interp_mesh( mdl, 0);
%           mdl_pts(i,:) is centre of element #i      
%   example: for mdl_pts= interp_mesh( mdl, 2);
%           mdl_pts(i,:,:) is 1 x [x,y,{z}] x n_points to interpolate
% 
% EXAMPLE:
%   mdl.nodes= [4,6,8,4,6,8;2,2,2,5,5,5]';
%   mdl.elems=[1,2,4;2,4,5;2,3,5;3,5,6];
%   mdl.type='fwd_model';mdl.name='jnk';
%   pp=interp_mesh(mdl,4);
%   show_fem(mdl); hold on ;
%      for i=1:size(pp,3); plot(pp(:,1,i),pp(:,2,i),'*'); end ;
%   hold off
% 

% (C) 2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if nargin<2; n_interp=0; end
[n_elems, n_dims_1]= size(mdl.elems);

% Get element nodes, and reshape
% need to be of size n_dims_1 x (n_elems*n_dims) for reshape
el_nodes= mdl.nodes(mdl.elems',:);
el_nodes= reshape(el_nodes, n_dims_1, []);

% Get interpolation matrix
interp= triangle_interpolation( n_interp, n_dims_1 );
l_interp = size(interp,1);

mdl_pts = interp*el_nodes;
mdl_pts = reshape(mdl_pts, l_interp, n_elems, n_dims_1-1);

mdl_pts = permute(mdl_pts, [2,3,1]);

% interpolate over a triangle with n_interp points
% generate a set of points to fairly cover the triangle
% dim_coarse is dimensions + 1 of coarse model
function interp= triangle_interpolation(n_interp, n_dims_1)
    interp= zeros(0,n_dims_1);

    if n_dims_1==3
       for i=0:n_interp
          for j=0:n_interp-i
             interp= [interp;i,j,n_interp-i-j];
          end
       end
    elseif n_dims_1==4
       for i=0:n_interp
          for j=0:n_interp-i
             for k=0:n_interp-i-j
                interp= [interp;i,j,k,n_interp-i-j-k];
             end
          end
       end
    else
       error('cant handle n_dims_1!=2');
    end

    interp= (interp + 1/n_dims_1 )/(n_interp+1);
