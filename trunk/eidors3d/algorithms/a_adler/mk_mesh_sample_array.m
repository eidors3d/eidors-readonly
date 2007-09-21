function [eptr, pt_coord]= mk_mesh_sample_array( varargin )
% MK_MESH_SAMPLE_ARRAY - create rectangular coordinate map into mesh
% [eptr, pt_coord]= mk_mesh_sample_array( mdl1, mdl2, ... , npoints);
%
% pt_coord - n_points x dims matrix of sample points
% eptr -     element number containing the pt_coord
% mdl1, mdl2 - fwd_model objects
%
% for example: [ep, pt] = mk_coarse_fine_mapping( m1, m2 )
% ep= [x   y   z]  pt= [m1, m2 ]
%      1   1   1         1  2 % pt in el#1 in m1, el#2 in m2
%      1   1   2         2  2 % pt in el#2 in m1, el#2 in m2
%      1   2   1         2  0 % pt in el#2 in m1, pt not in m2
%
% npoints - total number of points for rectangular mapping

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: mk_mesh_sample_array.m,v 1.3 2007-09-21 19:48:00 aadler Exp $

[pp, mdls] = proc_input (varargin);
xyz = interpxyz(pp.min_xyz , pp.max_xyz, pp.npoints);
tri_pts = mk_tri_pts( mdls{1}, xyz);
keyboard

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
    if any(csum<10)
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
   for mdl = mdls{:}'
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
