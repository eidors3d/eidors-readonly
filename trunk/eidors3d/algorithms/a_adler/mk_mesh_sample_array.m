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
% $Id: mk_mesh_sample_array.m,v 1.2 2007-09-21 14:18:36 aadler Exp $

[params, mdls] = proc_input (varargin);

keyboard

function [p, mdls] = proc_input (varargin);
   if isnumeric( varargin{1}{end} )
      p.npoints = varargin{1}{end};
      mdls = varargin{1}(1:end-1);
   else
      p.npoints = [];
      mdls = varargin{1}{:};
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

   if ~isempty(p.npoints)
      p.npoints = 64^p.dims;
   end
