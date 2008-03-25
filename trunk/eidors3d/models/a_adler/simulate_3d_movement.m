function [vh,vi,xyzr_pt]= simulate_3d_movement( n_sims, mdl_3d, rad_pr )
% SIMULATE_3D_MOVEMENT simulate rotational movement in 2D
% [vh,vi,xyzr_pt]= simulate_3d_movement( n_points, model, rad_pr )
%
% the target starts at (rad_pr(1),0) and rotates around 
%  clockwise
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
% OUTPUT:
%   vh - homogeneous measurements            M x 1
%   vi - target simulations                  M x n_points
%   xyzr_pt - x,y,z and radius of each point 3 x n_points

% $Id: simulate_3d_movement.m,v 1.5 2008-03-25 02:10:35 aadler Exp $

    if nargin <1
       n_sims = 200;
    end

    if nargin<2 || isempty(mdl_3d) % create our own fmdl
       mdl_3d= mk_common_model('n3r2',[16,2]); % NP's demo model
       mdl_3d= mdl_3d.fwd_model;
    end

    if nargin<3
       rad_pr= [2/3, 0.05, 0.1, 0.9];
    end


mdl_3d.solve=      'np_fwd_solve';
mdl_3d.system_mat= 'np_calc_system_mat';
mdl_3d.jacobian=   'np_calc_jacobian';
mdl_3d.misc.perm_sym=   '{n}';
n_elems= size(mdl_3d.elems,1);


% create homogeneous image + simulate data
sigma= ones( size(mdl_3d.elems,1) ,1);
img= eidors_obj('image', 'homogeneous image', ...
    'elem_data', sigma, ...
    'fwd_model', mdl_3d );
vh = fwd_solve( img);
% show_fem( homg_img);
eidors_msg('simulate_3d_movement: step #1: homogeneous simulation',2);

npx=128;
npy=128;
npz=64;
[x,y,z,radius,rp,z0,zt] = calc_point_grid(mdl_3d.nodes', npx, npy, npz, rad_pr);
clear pts;
for i=1:n_sims
    f_frac= (i-1)/n_sims;

    xp= radius * cos(f_frac*2*pi);
    yp= radius * sin(f_frac*2*pi);
    zp = z0 + (zt - z0) * f_frac;% object moves from bottow to top
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
eidors_msg('simulate_3d_movement: step #2: find points',2);

[eptr,vol]= img_mapper3a(mdl_3d.nodes', mdl_3d.elems',  ...
         x(pts_all), y(pts_all), z(pts_all));

target_conductivity= .2;

for i=1:n_sims
    obj_n= sparse( eptr(pts_idx{i}),1,1, n_elems, 1);
    img.elem_data= 1 + target_conductivity * obj_n./vol;
%   show_fem(img); view([-2,84]);pause;

    vi(i)= fwd_solve( img );% measurement
end
eidors_msg('simulate_3d_movement: step #3: target simulations',2);

vi= [vi(:).meas];
vh= [vh(:).meas];


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

function [x,y,z, radius, rp, zmin, zmax] = ...
         calc_point_grid(NODE, npx, npy, npz, rad_pr);

   xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
   xmean= mean([xmin,xmax]); xrange= xmax-xmin;

   ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
   ymean= mean([ymin,ymax]); yrange= ymax-ymin;

   zmin = min(NODE(3,:));    zmax = max(NODE(3,:));
   zmean= mean([zmin,zmax]); zrange= zmax-zmin;

   range= max([xrange, yrange,zrange]);
   [x y z]=ndgrid( ...
       linspace( xmean - range*0.5, xmean + range*0.5, npx ), ...
       linspace( ymean + range*0.5, ymean - range*0.5, npy ),...
       linspace( zmean - zrange*0.5, zmean + zrange*0.5, npz ));

   radius= rad_pr(1)*(xmax-xmin)/2;
   rp=     rad_pr(2)*(xmax-xmin)/2;
   zmin=   (rad_pr(3)-.5)*zrange + zmean;
   zmax=   (rad_pr(4)-.5)*zrange + zmean;
