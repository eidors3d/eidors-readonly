function [vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, mdl_3d, rad_pr,movefcn )
% SIMULATE_3D_MOVEMENT simulate rotational movement in 2D
% [vh,vi,xyzr_pt]= simulate_3d_movement( n_points, model, rad_pr, movefcn )
%
%   rad_pr = [path_radius, target_radius, zmin, zmax]
%      values are the fraction of the extent in each dimension
%      DEFAULT: [2/3, .05, .1, .9 ]
% 
%   n_points = number of points to simulate (default = 200)
%
%   model = fwd_model to simulate 
%         (default use internal, or if model= []);
%
%   movefcn = 1 (Default)  helical motion where the target starts
%     at (rad_pr(1),0) and rotates clockwise moving bottom to top.
%   movefcn = 2            radial movement in the vertical plane
%
%   movefcn = FUCN NAME or FUNC HANDLE
%      the function must accept the following parameters
%      [xp,yp,zp] = movefcn(f_frac, radius, z_bottom,z_top);
%
% OUTPUT:
%   vh - homogeneous measurements            M x 1
%   vi - target simulations                  M x n_points
%   xyzr_pt - x,y,z and radius of each point 3 x n_points

% $Id$

if nargin <1
   n_sims = 200;
end

if nargin<2 || isempty(mdl_3d) % create our own fmdl
   mdl_3d= mk_common_model('n3r2',[16,2]); % NP's demo model
   mdl_3d= mdl_3d.fwd_model;
end

if nargin<3 || isempty(rad_pr);
   rad_pr= [2/3, 0.05, 0.1, 0.9];
end

if nargin<4
   movefcn = 1;
end
if isnumeric(movefcn)
   if     movefcn==1
      movefcn = @helical_path;
   elseif movefcn==2
      movefcn = @radial_path;
   else
      error('value of movefcn not understood');
   end
else
   % assume movefcn is a function 
end


mdl_3d.solve=      'np_fwd_solve';
mdl_3d.system_mat= 'np_calc_system_mat';
mdl_3d.jacobian=   'np_calc_jacobian';
mdl_3d.misc.perm_sym=   '{n}';
n_elems= size(mdl_3d.elems,1);


eidors_msg('simulate_3d_movement: step #1: homogeneous simulation',2);
% create homogeneous image + simulate data
sigma= ones( size(mdl_3d.elems,1) ,1);
img= eidors_obj('image', 'homogeneous image', ...
    'elem_data', sigma, ...
    'fwd_model', mdl_3d );
vh = fwd_solve( img);

eidors_msg('simulate_3d_movement: step #2: find points',2);

if 0 % Old Style
   npx=128;
   npy=128;
   npz=64;
   [radius,rp,z0,zt,x,y,z] = calc_point_grid(mdl_3d.nodes', rad_pr, npx, npy, npz);

   clear pts;
   for i=1:n_sims
       f_frac= (i-1)/n_sims;

       % call function to simulate data
       [xp,yp,zp]= feval(movefcn, f_frac, radius, z0,zt);

       xyzr_pt(:,i)= [xp;-yp;zp;rp]; % -y because images and axes are reversed

       ff= find( (x(:)-xp).^2 + (y(:)-yp).^2 + (z(:)-zp).^2 <= rp^2 )';
       pts{i} = ff;
   end
   pts_all = unique( [pts{:}] );
   pts_all = pts_all(:);
   for i=1:n_sims
      [jnk,idx_i]= intersect( pts_all, pts{i});
      pts_idx{i}= idx_i;
   end

   [eptr,vol]= img_mapper3a(mdl_3d.nodes', mdl_3d.elems',  ...
            x(pts_all), y(pts_all), z(pts_all));
else
    mdl_pts = interp_mesh( mdl_3d, 4); % 45 per elem
    x= mdl_pts(:,1,:);
    y= mdl_pts(:,2,:);
    z= mdl_pts(:,3,:);
   [radius,rp,z0,zt] = calc_point_grid(mdl_3d.nodes', rad_pr);
end

target_conductivity= .1;

for i=1:n_sims
   if rem(i, max( floor(i/10), 10))==1; eidors_msg( ...
       'simulate_3d_movement: step #3 (%d of %d): target simulations', ...
       i, n_sims, 2); 
   end
   if 0
      obj_n= sparse( eptr(pts_idx{i}),1,1, n_elems, 1);
      img.elem_data= 1 + target_conductivity * obj_n./vol;
  %   show_fem(img); view([-2,84]);pause;
   else
      f_frac= (i-1)/n_sims;
      [xp,yp,zp]= feval(movefcn, f_frac, radius, z0,zt);
      xyzr_pt(:,i)= [xp;yp;zp;rp]; % -y because images and axes are reversed
      ff=  (x-xp).^2 + (y-yp).^2 + (z-zp).^2 <= rp^2;
      img.elem_data= 1 + target_conductivity * mean(ff,3);
   end

   vi(i)= fwd_solve( img );% measurement
end

vi= [vi(:).meas];
vh= [vh(:).meas];

%   movefcn = 1 (Default)  helical motion where the target starts
%     at (rad_pr(1),0) and rotates clockwise moving bottom to top.
% calculate x,y,z position of point, given f_frac of path
function [xp,yp,zp]= helical_path(f_frac, radius, z0,zt);
   xp= radius * cos(f_frac*2*pi);
   yp= radius * sin(f_frac*2*pi);
   % object moves from bottow to top
   zp = z0 + (zt - z0) * f_frac;

%   movefcn = 2            radial movement in the vertical plane
function [xp,yp,zp]= radial_path(f_frac, radius, z0,zt);
   rp= f_frac*radius; 
   cv= 2*pi*f_frac;
   xp= rp * cos(cv);
   yp= rp * sin(cv);
   zp = mean([zt,z0]);

% modified img_mapper3 from calc_slices.m
% this is like tsearch, but doesn't require
% delaunay triangularization. I also wrote it, so I like it :-)
function [EPTR, VOL] = img_mapper3a(NODE, ELEM, x,y,z );
   % calc and see if one point is in one element
   elems = size(ELEM,2);

   EPTR= zeros(prod(size(x)),1);
   VOL= zeros(elems,1);

   for j= 1: elems
       if rem(j,5000)==0; fprintf('mapping %d / %d \n',j,elems); end
       xyz= NODE(:,ELEM(:,j))';
       min_x= min(xyz(:,1)); max_x= max(xyz(:,1));
       min_y= min(xyz(:,2)); max_y= max(xyz(:,2));
       min_z= min(xyz(:,3)); max_z= max(xyz(:,3));

       % Simplex relative volume is det([v2-v1,v3-v1, ...])
       VOL(j)= abs(det(xyz'*[-1,1,0,0;-1,0,1,0;-1,0,0,1]'));

       xlmax= x<=max_x; if ~any(xlmax); continue; end
       xgmin= x>=min_x; if ~any(xgmin); continue; end
       ylmax= y<=max_y; if ~any(ylmax); continue; end
       ygmin= y>=min_y; if ~any(ygmin); continue; end
       zlmax= z<=max_z; if ~any(zlmax); continue; end
       zgmin= z>=min_z; if ~any(zgmin); continue; end
       % come up with a limited set of candidate points which
       % may be within the simplex
       endr=find( xlmax & xgmin & ylmax & ygmin & zlmax & zgmin);
       ll=  prod(size(endr));
       if isempty(ll);
          continue;
       end

       nn=  size(ELEM,1); %Simplex vertices
       vol=zeros(ll,nn);
       for i=1:nn
           i1= i; i2= rem(i,nn)+1; i3= rem(i+1,nn)+1;
           x1= xyz(i1,1)-x(endr); y1= xyz(i1,2)-y(endr); z1= xyz(i1,3)-z(endr);
           x2= xyz(i2,1)-x(endr); y2= xyz(i2,2)-y(endr); z2= xyz(i2,3)-z(endr);
           x3= xyz(i3,1)-x(endr); y3= xyz(i3,2)-y(endr); z3= xyz(i3,3)-z(endr);
           vol(:,i)= x1.*y2.*z3 - x1.*y3.*z2 - x2.*y1.*z3 + ...
               x3.*y1.*z2 + x2.*y3.*z1 - x3.*y2.*z1;
       end

       endr( sum(abs(vol),2) - VOL(j) >1e-8 )=[];
       EPTR(endr)= j;
   end %for j=1:ELEM

function [radius, rp, zmin, zmax,x,y,z] = ...
         calc_point_grid(NODE, rad_pr, npx, npy, npz);

   xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
   xmean= mean([xmin,xmax]); xrange= xmax-xmin;

   ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
   ymean= mean([ymin,ymax]); yrange= ymax-ymin;

   zmin = min(NODE(3,:));    zmax = max(NODE(3,:));
   zmean= mean([zmin,zmax]); zrange= zmax-zmin;

   radius= rad_pr(1)*(xmax-xmin)/2;
   rp=     rad_pr(2)*(xmax-xmin)/2;
   zmin=   (rad_pr(3)-.5)*zrange + zmean;
   zmax=   (rad_pr(4)-.5)*zrange + zmean;

   if nargout<=4; return; end
   range= max([xrange, yrange,zrange]);
   [x y z]=ndgrid( ...
       linspace( xmean - range*0.5, xmean + range*0.5, npx ), ...
       linspace( ymean + range*0.5, ymean - range*0.5, npy ),...
       linspace( zmean - zrange*0.5, zmean + zrange*0.5, npz ));
