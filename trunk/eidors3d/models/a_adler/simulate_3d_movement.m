function [vh,vi,xyr_pt]= simulate_3d_movement( n_sims, mdl_3d, rad_pr )
% SIMULATE_3D_MOVEMENT simulate rotational movement in 2D
% [vh,vi,xyr_pt]= simulate_3d_movement( n_points, model, rad_pr )
%
% the target starts at (rad_pr(1),0) and rotates around 
%  clockwise
% 
%   rad_pr = [path_radius, target_radius] = [2/3, .05] (default)
% 
%   n_points = number of points to simulate (default = 200)
%
%   model = fwd_model to simulate 
%         (default use internal, or if model= []);
%
% $Id: simulate_3d_movement.m,v 1.2 2008-03-24 14:55:01 aadler Exp $


mdl_3d.solve=      'np_fwd_solve';
mdl_3d.system_mat= 'np_calc_system_mat';
mdl_3d.jacobian=   'np_calc_jacobian';
mdl_3d.misc.perm_sym=   '{n}';
n_elems= size(mdl_3d.elems,1);


disp('STEP 1A: simultion 3D - homogeneous');
% create homogeneous image + simulate data
sigma= ones( size(mdl_3d.elems,1) ,1);
img= eidors_obj('image', 'homogeneous image', ...
    'elem_data', sigma, ...
    'fwd_model', mdl_3d );
vh = fwd_solve( img);
% show_fem( homg_img);

disp('STEP 1B: simultion 3D - inhomogeneous');
target_conductivity= .2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% np= calc_colours('npoints');
npx=64;
npy=64;
npz=64;
nn= mdl_3d.nodes;
[x,y,z]=ndgrid( linspace(min(nn(:,1)),max(nn(:,1)),npx), ...
                linspace(min(nn(:,2)),max(nn(:,2)),npy), ...
                linspace(min(nn(:,3)),max(nn(:,3)),npz));
% z= meshgrid(linspace(-.5,.5,np));
eptr= img_mapper3a(mdl_3d.nodes', mdl_3d.elems', npx, npy, npz);

eptr_n= sparse( eptr(:)+1,1,1, n_elems+1, 1);
eptr_n= full(eptr_n(2:end))+.01;%plus a small positive value to prevent providing 0

n_sims = 200;%frames number

rp= 1;%.05;%obj radius
radius=10;% 2/3;%obj trajactory radius (x-y plane)
z0 =5;% -.3;
zt =25;% .3;%obj moves along zaxis from z0 to zt
% rp= .05;%obj radius
% radius= 2/3;%obj trajactory radius (x-y plane)
% z0 = -.3;
% zt = .3;%obj moves along zaxis from z0 to zt

for i=1:n_sims
    f_frac= (i-1)/n_sims;
    if rem(i,20)==0; fprintf('simulating %d / %d \n',i,n_sims); end

    xp= radius * cos(f_frac*2*pi);
    yp= radius * sin(f_frac*2*pi);
    zp = z0 + (zt - z0) * f_frac;% object moves from bottow to top
    xyzr_pt(:,i)= [xp;-yp;zp;rp]; % -y because images and axes are reversed

    ff= find( (x(:)-xp).^2 + (y(:)-yp).^2 + (z(:)-zp).^2 <= rp^2 )';%find points in obj
    obj_n= sparse( eptr(ff)+1,1,1, n_elems+1, 1);
    obj_n= full(obj_n(2:end));

    img.elem_data= 1 + target_conductivity * (obj_n./eptr_n);%inhom img

    vi(i)= fwd_solve( img );% measurement
end

vi= [vi(:).meas];
vh= [vh(:).meas];


% modified img_mapper3 from calc_slices.m
% find out indices of points that locate in elements
function EPTR= img_mapper3a(NODE, ELEM, npx, npy,npz )

   [x,y,z] = calc_point_grid(NODE, ELEM, npx, npy, npx);

   EPTR=zeros(npy,npx,npz);

   % calc and see if one point is in one element
   elems = size(ELEM,2);
   for j= 1: elems
       if rem(j,1000)==0; fprintf('mapping %d / %d \n',j,elems); end
       xyz= NODE(:,ELEM(:,j))';
       min_z= min(xyz(:,3)); max_z= max(xyz(:,3));
       %     if (min_z>0 | max_z<0)
       %         continue;
       %     end
       min_x= min(xyz(:,1)); max_x= max(xyz(:,1));
       min_y= min(xyz(:,2)); max_y= max(xyz(:,2));

       % Simplex volume is det([v2-v1,v3-v1, ...])
       VOL= abs(det(xyz'*[-1,1,0,0;-1,0,1,0;-1,0,0,1]'));

       % come up with a limited set of candidate points which
       % may be within the simplex
       endr=find( y(:)<=max_y & y(:)>=min_y ...
                & x(:)<=max_x & x(:)>=min_x ...
                & z(:)<=max_z & z(:)>=min_z);

       nn=  size(ELEM,1); %Simplex vertices
       ll=  length(endr);
       vol=zeros(ll,nn);
       for i=1:nn
           i1= i; i2= rem(i,nn)+1; i3= rem(i+1,nn)+1;
           x1= xyz(i1,1)-x(endr); y1= xyz(i1,2)-y(endr); z1= xyz(i1,3)-z(endr);
           x2= xyz(i2,1)-x(endr); y2= xyz(i2,2)-y(endr); z2= xyz(i2,3)-z(endr);
           x3= xyz(i3,1)-x(endr); y3= xyz(i3,2)-y(endr); z3= xyz(i3,3)-z(endr);
           vol(:,i)= x1.*y2.*z3 - x1.*y3.*z2 - x2.*y1.*z3 + ...
               x3.*y1.*z2 + x2.*y3.*z1 - x3.*y2.*z1;
       end

       endr( sum(abs(vol),2) - VOL >1e-8 )=[];
       EPTR(endr)= j;
   end %for j=1:ELEM

function [x,y,z] = calc_point_grid(NODE, ELEM, npx, npy, npx);
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
