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
%   n_interp may be specified as:
%      fwd_model.interp_mesh.n_interp (This overrides the above)
%
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

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test;return; end

if nargin<2; n_interp=0; end
try n_interp = mdl.interp_mesh.n_interp; end % Override if provided


% cashing
   
   c_obj = {mdl.elems, mdl.nodes, n_interp};
   mdl_pts = eidors_obj('get-cache', c_obj, 'interpolation');
   if ~isempty(mdl_pts)
       return
   end


% Get element nodes, and reshape
% need to be of size n_dims_1 x (n_elems*n_dims) for reshape
el_nodes= mdl.nodes(mdl.elems',:);
el_nodes= reshape(el_nodes, elem_dim(mdl)+1, []);

% Get interpolation matrix
interp= triangle_interpolation( n_interp, elem_dim(mdl) );
l_interp = size(interp,1);

mdl_pts = interp*el_nodes;
mdl_pts = reshape(mdl_pts, l_interp, num_elems(mdl), mdl_dim(mdl));

mdl_pts = permute(mdl_pts, [2,3,1]);

% caching
eidors_cache('boost_priority', -2); % low priority
c_obj = {mdl.elems, mdl.nodes, n_interp};
eidors_obj('set-cache', c_obj, 'interpolation', mdl_pts);
eidors_cache('boost_priority', +2); % restore priority


% interpolate over a triangle with n_interp points
% generate a set of points to fairly cover the triangle
% dim_coarse is dimensions + 1 of coarse model
function interp= triangle_interpolation(n_interp, el_dim)
    interp= zeros(0,el_dim+1);

    if el_dim==2
       for i=0:n_interp
          for j=0:n_interp-i
             interp= [interp;i,j,n_interp-i-j];
          end
       end
    elseif el_dim==3
       for i=0:n_interp
          for j=0:n_interp-i
             for k=0:n_interp-i-j
                interp= [interp;i,j,k,n_interp-i-j-k];
             end
          end
       end
    else
       error('Element is not 2D (triangle) or 3D (tetrahedron)');
    end

    interp= (interp + 1/(el_dim+1) )/(n_interp+1);


function do_unit_test
    mdl.nodes= 3*[4,6,8,4,6,8;2,2,2,5,5,5]';
    mdl.elems=[1,2,4;2,4,5;2,3,5;3,5,6];
    mdl.type='fwd_model';mdl.name='jnk';
    pp=interp_mesh(mdl,0);
    unit_test_cmp('2D/2D (#1): ',pp,[14 9;16 12; 20 9;22 12]);

    pp=interp_mesh(mdl,3);
    unit_test_cmp('2D/2D (#2): ',pp(:,:,6),[14 9;16 12; 20 9;22 12]);

    mdl.nodes = [mdl.nodes, 0*mdl.nodes(:,1)+3];
    pp=interp_mesh(mdl,0);
    unit_test_cmp('2D/3D (#1): ',pp,[14 9 3;16 12 3; 20 9 3;22 12 3]);

    mdl = mk_common_model('n3r2',16); mdl= mdl.fwd_model;
    pp=interp_mesh(mdl,0);
    unit_test_cmp('3D/3D (#1): ',pp(1,:),[0.920196320100808   0.048772580504032   0.5],1e-14);

    pp=interp_mesh(mdl,4);
    unit_test_cmp('3D/3D (#2a):',pp(1,:,21),[0.920196320100808   0.048772580504032   0.5],1e-14);

