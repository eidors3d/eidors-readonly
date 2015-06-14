function quiv = show_current( img, vv )
% SHOW_CURRENT: show current or other quantity defined
%  on nodes onto the image
%
% show_current( img, volt )
%   img -> img object 
%   volt-> voltage on nodes (if not specified, img is solved
%      via fwd_solve, or the value on node_data is used)
%
% Without specified values, show_current will
%   create one current vector for each element
% The points are specified as either
%   img.fwd_model.mdl_slice_mapper.npx   - number of points in horizontal direction
%   img.fwd_model.mdl_slice_mapper.npy   - number of points in vertical 
%    or
%   img.fwd_model.mdl_slice_mapper.x_pts - vector of points in horizontal direction
%   img.fwd_model.mdl_slice_mapper.y_pts - vector of points in vertical
%
% For 3D models, the slice may be specified as (see mdl_slice_mapper for information)
%   img.fwd_model.mdl_slice_mapper.level 
% 
% If an output is specified, then no image is draws
% 
% Examples:
%   img = mk_image( mk_common_model('b2c2',8));
%   show_current(img); 
% OR
%   img = mk_image( mk_common_model('b2c2',8));
%   img.fwd_model.mdl_slice_mapper.npx  = 64;
%   img.fwd_model.mdl_slice_mapper.npy  = 32;
%   show_current(img); 
% OR
%   img = mk_image( mk_common_model('b2c2',8));
%   q = show_current(img); 
%   quiver(q.xp,q.yp, q.xc,q.yc, 2,'b','LineWidth',2);


% (C) 2010 Andy Adler. License: GPL version 2 or version 3
% $Id$

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

dims = size(img.fwd_model.nodes,2);

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
  if dims ==2
     VE  = ([oo, nod])\fix_dim(vv(idx));
  else
     VE  = ([oo, nod])\vv(idx);
  end
  elemcur(i+1,:) = img.elem_data(i,1)*VE(2:end)';
%  elemcur(i+1,:) = (reshape(img.elem_data(i,1,:,:),dims,dims)*VE(2:end))';

end

if isfield(img.fwd_model, 'mdl_slice_mapper');
   if dims == 2;
      img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   end
   elem_ptr = mdl_slice_mapper( img.fwd_model, 'elem' );
   szep = size(elem_ptr);

   [xp,yp] = grid_the_space( img.fwd_model);

   xc = reshape( elemcur(elem_ptr+1,1), szep);
   yc = reshape( elemcur(elem_ptr+1,2), szep);
   if dims==3
      zc = reshape( elemcur(elem_ptr+1,3), szep);
   end
else 
   pp = interp_mesh( img.fwd_model);
   xp = pp(:,1);
   yp= pp(:,2);
   xc = elemcur(2:end,1);
   yc = elemcur(2:end,2);
   if dims==3
      zc = elemcur(2:end,3);
   end
end

quiv.xp = xp; 
quiv.yp = yp; 
quiv.xc = xc; 
quiv.yc = yc; 
if dims==3
   quiv.zc = zc; 
end
if nargout == 0
   quiver(xp,yp,xc,yc,2,'k');
end

function vv = fix_dim(vv)
    if size(vv,1) == 1
        vv = vv';
    end

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
  if size(xspace,2) == 1
      [x,y]=meshgrid( xspace, yspace );
  else
      x= xspace;
      y= yspace;
  end

function do_unit_test
   fmdl.nodes = [0,0;0,1;1,0;1,1];
   fmdl.elems = [1,2,3;2,3,4];
   fmdl.electrode(1).nodes = [1,2]; fmdl.electrode(1).z_contact = 0.01;
   fmdl.electrode(2).nodes = [3,4]; fmdl.electrode(2).z_contact = 0.01;
   fmdl.gnd_node = 1;
   fmdl.stimulation(1).stim_pattern = [1;-1];
   fmdl.stimulation(1).meas_pattern = [1,-1];
   fmdl.solve = @fwd_solve_1st_order;
   fmdl.system_mat = @system_mat_1st_order;
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
