function fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, spacing);
% DM_2D_CIRC_PT_ELECS: Create circle mesh (or radius 1) refined with electrodes
% fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, spacing);
%     at points on the electrodes
% fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, ...
%               [base_spacing, refine_ratio, gradient]);
% elec_pts = cell fcn of N x [x,y,{z}] for each electrode
%    normally only two points are specified (at start and end of electrode)
%    eg elec_pts{1} = [-.1,1;.1,1];
% pfix = any fixed points to provide to the model (default = [])
% param(1) = base_spacing - edge length away from refined nodes (eg 0.1)
% param(2) = refine_ratio - relative refinement near points (eg. 10)
% param(3) = gradient     - transition slope of refinement (eg 0.1)
%
% Example: Three electrodes 2 on boundary and one internal
%  elec_pts = {[1,0],[0,1;sin(0.2),cos(0.2)],[0.5,0.5]};
%  fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.15,10,0.05] );
%
% Example:
%  n_elecs= 14; elec_width= 0.1; hw= elec_width/2;
%  th = linspace(0,2*pi,n_elecs+1); th(end)=[];
%  for i=1:n_elecs;
%     ti = th(i) + [hw;-hw];
%     elec_pts{i} = [sin(ti),cos(ti)];
%  end
%  fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );

%
% See also: dm_2d_pt_elecs

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

cache_obj = {elec_pts, spacing};
fmdl= eidors_obj('get-cache',cache_obj, 'dm_2d_circ_pt_elecs');
if ~isempty(fmdl); return; end % Object is cached

params.base_spacing = spacing(1);
params.refine_ratio = spacing(2);
params.gradient     = spacing(3);

bbox= [-1,-1;1,1];
fmdl= dm_2d_pt_elecs( elec_pts, [], params, @circle, [-1,-1;1,1] );

fmdl.name = sprintf('dm_2d_circ_pt_elec');
eidors_obj('set-cache',cache_obj, 'dm_2d_circ_pt_elecs', fmdl);

function d= circle(p,params);
  d = sqrt(sum(p.^2,2)) - 1; 
