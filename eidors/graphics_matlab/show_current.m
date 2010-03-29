function show_current( img, vv )
% SHOW_CURRENT: show current or other quantity defined
%  on nodes onto the image
%
% show_current( img, volt )
%   img -> img object 
%   volt-> voltage on nodes (if not specified, img is solved via fwd_solve)
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
   img.fwd_solve.get_all_meas = 1;
   vh = fwd_solve(img);
   vv = vh.volt(:,1);
end 

elem_ptr = mdl_slice_mapper( img.fwd_model, 'elem' );
szep = size(elem_ptr);

Nel = size(img.fwd_model.elems,1);
elemcur = zeros(Nel+1,dims);
del = [1,-1,0;0,1,-1;-1,0,1];
del2 = del'*del;
for i=1:Nel
  idx = img.fwd_model.elems(i,:);
  nod = img.fwd_model.nodes(idx,:);
  elemcur(i+1,:) = vv(idx)'*del2*nod;
end

[xp,yp] = grid_the_space( img.fwd_model);

xc = reshape( elemcur(elem_ptr+1,1), szep);
yc = reshape( elemcur(elem_ptr+1,2), szep);

quiver(xp,yp,xc,yc,2);

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
