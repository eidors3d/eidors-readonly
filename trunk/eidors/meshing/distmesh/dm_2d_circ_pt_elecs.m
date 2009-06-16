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


params.base_spacing = spacing(1);
params.refine_ratio = spacing(2);
params.gradient     = spacing(3);
epfix = vertcat(elec_pts{:});
params.refine_pts   = epfix;

[p,t] = distmeshnd(@circle,@dm_refine_points,params.base_spacing/2,  ...
           [-1,-1;1,1],[pfix,epfix],params);

bdy = find_boundary(t);
ubn = unique(bdy(:));

for i=1:length(elec_pts);
   electrode(i).nodes     = get_elec_nodes(elec_pts{i}, p, ubn);
   electrode(i).z_contact = 0.01; % nominal
end


fmdl.nodes= p;
fmdl.elems= t;
fmdl.boundary = bdy;
fmdl.type = 'fwd_model';
fmdl.name = sprintf('dm_2d_circ_pt_elec');
fmdl.electrode = electrode;
fmdl.gnd_node=           1;
fmdl.np_fwd_solve.perm_sym =     '{n}';


function d= circle(p,params);
  d = sqrt(sum(p.^2,2)) - 1; 

% Find the nodes associated with a given electrode
%   electrode(i).nodes = get_elec_nodes(elec_pts{i}, p, ubn);
% electrode points are no further from each end point
%   than the elecs are from each other
function nodes= get_elec_nodes( epts, p, ubn);
   tol = 1e-6;
   nodes = ubn;
   bp    = p(ubn,:); % Unique points on the boundary
   lep = size(epts,1);
   oep = ones(lep,1);
   onp = ones(size(ubn,1),1);
   for i=1:lep
      el = epts(i,:);
      % Distance to other electrodes
      de = oep*el - epts;
      de = sqrt( sum( de.^2,2 ));
      mde= max(de);
      % Distance to other points
      dp = onp*el - bp;
      dp = sqrt( sum( dp.^2,2 ));
      % Zero outside points
      nodes( dp > mde+tol ) = 0;
   end
   nodes = nodes( nodes>0 );
