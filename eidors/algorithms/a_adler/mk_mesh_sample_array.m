function [eptr, pt_xyz]= mk_mesh_sample_array( varargin )
% MK_MESH_SAMPLE_ARRAY - create rectangular coordinate map into mesh
% [eptr, pt_xyz]= mk_mesh_sample_array( mdl1, mdl2, ... , npoints);
%
% pt_xyz - n_points x dims matrix of sample points
% eptr -     element number containing the pt_coord
% mdl1, mdl2 - fwd_model objects
%
% for example: [ep, pt] = mk_coarse_fine_mapping( m1, m2 )
% ep= [x   y   z]  pt= [m1, m2 ]
%      1   1   1        1   2  % pt in el#1 in m1, el#2 in m2
%      1   1   2        2   2  % pt in el#2 in m1, el#2 in m2
%      1   2   1        2  NaN % pt in el#2 in m1, pt not in m2
%
% npoints - total number of points for rectangular mapping

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

[pp, mdls] = proc_input (varargin);
pt_xyz = interpxyz(pp.min_xyz , pp.max_xyz, pp.npoints);
for mdl= mdls(:)'
   pt_xyz= [pt_xyz; interpmdl( mdl{1} )];
end


eptr= [];
for mdl= mdls(:)'
   tri_pts = mk_tri_pts( mdl{1}, pt_xyz);
   eptr = [eptr, tri_pts];
end

% Add four interpolation points in each elem
% 1- n1+n2+n3/3 = mid
% 2- n1+mid/2 = n1/2 + (n1+n2+n3)/6 = (4*n1 + n2 + n3)/6;
function xyz = interpmdl( fmdl );
    to3d = speye(size(fmdl.nodes,2), 3);
    nodes = fmdl.nodes*to3d;
    elems = fmdl.elems;

    xelems= reshape( nodes(elems,1), size(elems) );
    yelems= reshape( nodes(elems,2), size(elems) );
    zelems= reshape( nodes(elems,3), size(elems) );

    if size(elems,2)==4 % pyramids
       iv= [2;2;2;2]/8;
       xyz=[      [xelems*iv, yelems*iv, zelems*iv]];

       iv= [5;1;1;1]/8;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];

       iv= [1;5;1;1]/8;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];

       iv= [1;1;5;1]/8;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];

       iv= [1;1;1;5]/8;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];
    elseif size(elems,2)==3 % triangles
       iv= [2;2;2]/6;
       xyz=[      [xelems*iv, yelems*iv, zelems*iv]];

       iv= [4;1;1]/6;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];

       iv= [1;4;1]/6;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];

       iv= [1;1;4]/6;
       xyz=[xyz;  [xelems*iv, yelems*iv, zelems*iv]];
    else
       error('strange size elems'); 
    end

function xyz = interpxyz( xyzmin, xyzmax, n_interp);
    xyz_delta= xyzmax - xyzmin;
    % Allocate points to each dir nx,ny,nz so nx*ny*nz~=n_interp,
    % and (nx,ny,nz) = k*xyzdelta
    neq0 = xyz_delta~=0;
    k= (prod(xyz_delta(neq0))/n_interp)^(1/sum(neq0));
    n_xyz= [1,1,1]; n_xyz(neq0)= ceil(xyz_delta(neq0)/k);

    xspace = linspace(xyzmin(1), xyzmax(1), n_xyz(1) );
    yspace = linspace(xyzmin(2), xyzmax(2), n_xyz(2) );
    zspace = linspace(xyzmin(3), xyzmax(3), n_xyz(3) );
    [xx3,yy3,zz3] = ndgrid( xspace, yspace, zspace );
    xyz= [xx3(:), yy3(:), zz3(:)];

% calculate mapping of points xyz into triangle
% fmdl is fwd_mode, xyz is points [x(:), y(:), z(:)]
function tri_pts = mk_tri_pts( fmdl, xyz);
    to3d = speye(size(fmdl.nodes,2), 3);
    nodes = fmdl.nodes*to3d;
    elems = fmdl.elems;
    tri_pts= tsearchn(nodes, elems, xyz*to3d');

    ff= find(~isnan(tri_pts));
    n_ff= size(elems, 1);
    csum = sparse(tri_pts(ff),1,1,n_ff, 1);
    if any(csum<4)
       warning(sprintf( ...
         'mesh density is too low: min=%d', full(min(csum)) ))
    end

function [p, mdls] = proc_input (varargin);
   % BLOODY STUPID IDIOTS AT MATLAB CANT KEEP COMPATIBLE
   varargin = varargin{1};

   if isnumeric( varargin{end} )
      p.npoints = varargin{end};
      mdls = varargin(1:end-1);
   else
      p.npoints = [];
      mdls = varargin(:);
   end

   p.dims= 1;
   p.min_xyz =  [inf,inf,inf];
   p.max_xyz = -[inf,inf,inf];
   for mdlc = mdls(:)'
      mdl= mdlc{1}; 
      try if ~strcmp( mdl.type , 'fwd_model');
         error('model is not a fwd_model');
      end; catch
         error('model is not a fwd_model');
      end

      nodes= mdl.nodes;

      max_n= [0,0,0];
      dim = size(nodes,2);
      if dim>p.dims; p.dims= dim; end

      max_n(1:dim) = max(nodes);
      idx = p.max_xyz < max_n;
      p.max_xyz(idx) = max_n(idx);

      min_n= [0,0,0];
      min_n(1:dim) = min(nodes);
      idx = p.min_xyz > min_n;
      p.min_xyz(idx) = min_n(idx);
   end

   if isempty(p.npoints)
      p.npoints = 64^p.dims;
   end
