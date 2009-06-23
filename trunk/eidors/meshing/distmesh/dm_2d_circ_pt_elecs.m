function fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, spacing);
% DM_2D_CIRC_PT_ELECS: Create circle mesh (or radius 1) refined
%     at points on the electrodes
% fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, ...
%               [base_spacing, refine_ratio, gradient]);
% elec_pts = cell fcn of N x [x,y,{z}] for each electrode
%    normally only two points are specified (at start and end of electrode)
%    eg elec_pts{1} = [-.1,1;.1,1];
% pfix = any fixed points to provide to the model (default = [])
% base_spacing - edge length away from refined nodes (eg 0.1)
% refine_ratio - relative refinement near points (eg. 10)
% gradient     - transition slope of refinement (eg 0.1)
%
% Example: = dm_2d_circ_pt_elecs( elec_pts, [], [0.15,10,0.05] );

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id$

bbox= [-1,-1;1,1];
fmdl= dm_2d_pt_elecs( elec_pts, [], [0.15,10,0.05], @circle, [-1,-1;1,1] );

fmdl.name = sprintf('dm_2d_circ_pt_elec');

function d= circle(p,params);
  d = sqrt(sum(p.^2,2)) - 1; 
