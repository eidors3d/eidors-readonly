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


if nargin==1;
  e_curr = calc_elem_current( img );
else 
  e_curr = calc_elem_current( img, vv);
end

dims = size(img.fwd_model.nodes,2);

elemcur = [zeros(1,dims); e_curr ];
try
   if strcmp(img.calc_colours.component, 'real')
       elemcur = real(elemcur);
   end
   if strcmp(img.calc_colours.component, 'imag')
       elemcur = imag(elemcur);
   end
end  % if no calc_colours.component field, do nothing

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
   fmdl.type = 'fwd_model';
   fmdl.normalize_measurements= 0;
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
