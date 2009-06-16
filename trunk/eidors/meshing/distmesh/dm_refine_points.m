function h = dm_refine_points( pts, params );
% DM_REFINE_POINTS: refine distmesh volume at point locations:
% h= dm_refine_points( pts, params );
%   pts is an array NxDims of node positions
% params.refine_pts   - points at which to refine mesh (NxNdims)
% params.base_spacing - edge length away from refined nodes (eg 0.1)
% params.refine_ratio - relative refinement near points (eg. 10)
% params.gradient     - transition slope of refinement (eg 0.1)

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

ep     = params.refine_pts;
maxsize= params.base_spacing;
minsize= params.base_spacing / params.refine_ratio;
grad   = params.gradient;

op =  ones(size(pts,1),1);
for i=1:size(ep,1);
   de = sqrt( sum( (pts - op*ep(i,:)).^2 ,2) );
   h_i = min(minsize+grad*de,maxsize);
   if i==1;   h=h_i;
   else       h = min(h,h_i);
   end
end

