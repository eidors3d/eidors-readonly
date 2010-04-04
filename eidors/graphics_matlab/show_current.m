function show_current( img, vv )
% SHOW_CURRENT: show current or other quantity defined
%  on nodes onto the image
%
% show_current( img, volt )
%   img -> img object 
%   volt-> voltage on nodes (if not specified, img is solved
%      via fwd_solve, or the value on node_data is used)
%
% The points are specified as either
%   img.fwd_model.mdl_slice_mapper.npx   - number of points in horizontal direction
%   img.fwd_model.mdl_slice_mapper.npy   - number of points in vertical 
%    or
%   img.fwd_model.mdl_slice_mapper.x_pts - vector of points in horizontal direction
%   img.fwd_model.mdl_slice_mapper.y_pts - vector of points in vertical
%
% For 3D models, the slice may be specified as (see calc_slices for information)
%   img.fwd_model.mdl_slice_mapper.level 

% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

dims = size(img.fwd_model.nodes,2);
if dims == 2;
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
end

if nargin==1; % We need to calculate
   if isfield(img,'elem_data')
      img.fwd_solve.get_all_meas = 1;
      vh = fwd_solve(img);
      vv = vh.volt(:,1);
   elseif isfield(img,'node_data');
      vv = img.node_data(:,1);
      error('show_current: cannot interpolate conductivity onto elements (yet)');
   else
      error('show_current: one parameter provided, and cannot solve for node voltages');
   end
end 

elem_ptr = mdl_slice_mapper( img.fwd_model, 'elem' );
szep = size(elem_ptr);

Nel = size(img.fwd_model.elems,1);
elemcur = zeros(Nel+1,dims);
% Calc field as I = sigma E
%V1 = V0 + Ex*x1 + Ey*y1   [ 1 x1 y1 ] [ V0 ]
%V2 = V0 + Ex*x2 + Ey*y2 = [ 1 x2 y2 ]*[ Ex ]
%V3 = V0 + Ex*x3 + Ey*y    [ 1 x3 y3 ] [ Ey ]
oo = ones(dims+1,1);
for i=1:Nel
  idx = img.fwd_model.elems(i,:);
  nod = img.fwd_model.nodes(idx,:);
  VE  = ([oo, nod])\vv(idx);
   
  elemcur(i+1,:) = img.elem_data(i)*VE(2:end)';
end

[xp,yp] = grid_the_space( img.fwd_model);

xc = reshape( elemcur(elem_ptr+1,1), szep);
yc = reshape( elemcur(elem_ptr+1,2), szep);

quiver(xp,yp,xc,yc,2,'k');

function  [x,y] = grid_the_space( fmdl);

  xspace = []; yspace = [];
  try 
     xspace = fmdl.mdl_slice_mapper.x_pts;
     yspace = fmdl.mdl_slice_mapper.y_pts;
  end

  if isempty(xspace)
     npx  = fmdl.mdl_slice_mapper.npx;
     npy  = fmdl.mdl_slice_mapper.npy;

     xmin = min(fmdl.nodes(:,1));    xmax = max(fmdl.nodes(:,1));
     xmean= mean([xmin,xmax]); xrange= xmax-xmin;

     ymin = min(fmdl.nodes(:,2));    ymax = max(fmdl.nodes(:,2));
     ymean= mean([ymin,ymax]); yrange= ymax-ymin;

     range= max([xrange, yrange]);
     xspace = linspace( xmean - range*0.5, xmean + range*0.5, npx );
     yspace = linspace( ymean + range*0.5, ymean - range*0.5, npy );
  end

  [x,y]=meshgrid( xspace, yspace );

function do_unit_test
   fmdl.nodes = [0,0;0,1;1,0;1,1];
   fmdl.elems = [1,2,3;2,3,4];
   fmdl.electrode(1).nodes = [1,2]; fmdl.electrode(1).z_contact = 0.01;
   fmdl.electrode(2).nodes = [3,4]; fmdl.electrode(2).z_contact = 0.01;
   fmdl.gnd_node = 1;
   fmdl.stimulation(1).stim_pattern = [1;-1];
   fmdl.stimulation(1).meas_pattern = [1,-1];
   fmdl.solve = @aa_fwd_solve;
   fmdl.system_mat = @aa_calc_system_mat;
   fmdl.type = 'fwd_model'
   img = mk_image(fmdl,[1,1]); 
   img.fwd_solve.get_all_meas = 1;

   img.fwd_model.mdl_slice_mapper.npx = 6;
   img.fwd_model.mdl_slice_mapper.npy = 6;
   show_current(img);

   imdl= mk_common_model('d2d2c',8);
   img = calc_jacobian_bkgnd( imdl );
   img.fwd_model.mdl_slice_mapper.npx = 64;
   img.fwd_model.mdl_slice_mapper.npy = 64;
   show_current(img);

   img.fwd_model.mdl_slice_mapper.x_pts = linspace(0,1,62);
   img.fwd_model.mdl_slice_mapper.y_pts = linspace(0,1,56);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];

   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve(img);
   show_current(img, vh.volt(:,1));
