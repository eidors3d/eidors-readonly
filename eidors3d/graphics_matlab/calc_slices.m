function rimg = calc_slices( img, levels );
% calc_slices (img, levels, clim  ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or array of I images
% levels = Matrix [Lx3] of L image levels
%          each row of the matrix specifies the intercepts
%          of the slice on the x, y, z axis. To specify a z=2 plane
%          parallel to the x,y: use levels= [inf,inf,2]
%
% if levels is scalar, then make levels equispaced horizontal
%          cuts through the object
%
% rimg= np x np x I x L where np is 128 by default
% np can be adjusted by calc_colours('npoints')

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_slices.m,v 1.21 2008-02-29 22:56:27 aadler Exp $

np= calc_colours('npoints');

fwd_model= img(1).fwd_model; % Assume all fwd_models are same
dims= size(fwd_model.nodes,2);
if nargin<=1;
   levels= [];
end

if isempty(levels) && dims==2
   levels= [Inf,Inf,0];
end

if size(levels)== [1,1]
   zmax= max(fwd_model.nodes(:,3));
   zmin= min(fwd_model.nodes(:,3));
   levels = linspace(zmin,zmax, levels+2);
   levels = levels(2:end-1)'*[Inf,Inf,1];
end

num_levs= size(levels,1);
if isfield(img,'elem_data')
   elem_data= [img.elem_data];
   if size(elem_data,1)==1; elem_data=elem_data';end
   n_images= size(elem_data,2);
   rimg=zeros(np,np,n_images,num_levs);

   for lev_no = 1:num_levs
      level= levels( lev_no, 1:3 );

      rimg(:,:,:,lev_no) = calc_image_elems( elem_data, level, fwd_model, np);
   end
elseif isfield(img,'node_data')
else
   error('img does not have a data field');
end


% Calculate an image by mapping it onto the elem_ptr matrix
function rimg= calc_image_elems( elem_data, level, fwd_model, np)

% elem_ptr_table also depends on the number of mapped points
fwd_model.calc_slices.mapping_npoints=np;

% Get elem_ptr from cache, if available
% EPtable is cell array of elem_ptrs for different levels
EPtable = eidors_obj('get-cache', fwd_model, 'elem_ptr_table');
level_hash= eidors_var_id( level );

if isfield(EPtable, level_hash);
   elem_ptr= getfield(EPtable, level_hash);
else
   [NODE, ELEM] = level_model( fwd_model, level );
   if size(NODE,1) ==2 %2D
      elem_ptr= img_mapper2( NODE, ELEM, np, np);
   else
      elem_ptr= img_mapper3( NODE, ELEM, np, np);
   end

   EPtable = setfield(EPtable, level_hash, elem_ptr);
   eidors_obj('set-cache', fwd_model, 'elem_ptr_table', EPtable);
   eidors_msg('show_slices: setting cached value', 3);
end


backgnd= NaN;
n_images= size(elem_data,2);
rval= [backgnd*ones(1,n_images); elem_data];
rimg= reshape( rval(elem_ptr+1,:), np,np, n_images );


% Search through each element and find the points which
% are in that element
function EPTR= img_mapper2(NODE, ELEM, npx, npy );
  xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
  xmean= mean([xmin,xmax]); xrange= xmax-xmin;

  ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  range= max([xrange, yrange]);
  [x y]=meshgrid( ...
      linspace( xmean - range*0.5, xmean + range*0.5, npx ), ...
      linspace( ymean + range*0.5, ymean - range*0.5, npy ) );
  v_yx= [-y(:) x(:)];
  turn= [0 -1 1;1 0 -1;-1 1 0];
  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A
  for j= 1: size(ELEM,2)
    % calculate area of three subtrianges to each candidate point.
    xy= NODE(:,ELEM(:,j))';
    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
             & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) );
    % a is determinant of matrix [i,j,k, xy]
    a= xy([2;3;1],1).*xy([3;1;2],2)- xy([3;1;2],1).*xy([2;3;1],2);

    aa= sum(abs(ones(length(endr),1)*a'+ ...
                v_yx(endr,:)*xy'*turn)');
    endr( abs( (abs(sum(a))-aa) ./ sum(a)) >1e-8)=[];
    EPTR(endr)= j;
  end %for j=1:ELEM

% 2D mapper of points to elements. First, we assume that
% The vertex geometry (NODE) has been rotated and translated
% so that the imaging plane is on the z-axis. Then we iterate
% through elements to find the containing each pixel
function EPTR= img_mapper2a(NODE, ELEM, npx, npy );
  xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
  xmean= mean([xmin,xmax]); xrange= xmax-xmin;

  ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  range= max([xrange, yrange]);
  [x y]=meshgrid( ...
      linspace( xmean - range*0.5, xmean + range*0.5, npx ), ...
      linspace( ymean + range*0.5, ymean - range*0.5, npy ) );

  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A
  for j= 1: size(ELEM,2)
    xyz= NODE(:,ELEM(:,j))';
    min_x= min(xyz(:,1)); max_x= max(xyz(:,1));
    min_y= min(xyz(:,2)); max_y= max(xyz(:,2));

    % Simplex volume is det([v2-v1,v3-v1, ...])
    VOL= abs(det(xyz'*[-1,1,0;-1,0,1]'));

    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max_y & y(:)>=min_y ...
             & x(:)<=max_x & x(:)>=min_x );

    nn=  size(ELEM,1); %Simplex vertices
    ll=  length(endr);
    vol=zeros(ll,nn);
    for i=1:nn
       i1= i; i2= rem(i,nn)+1;
       x1= xyz(i1,1) - x(endr);
       y1= xyz(i1,2) - y(endr);
       x2= xyz(i2,1) - x(endr);
       y2= xyz(i2,2) - y(endr);
       vol(:,i)= x1.*y2 - x2.*y1;  % determinant
    end

    endr( sum(abs(vol),2) - VOL >1e-8 )=[];
    EPTR(endr)= j;
  end %for j=1:ELEM


% 3D mapper of points to elements. First, we assume that
% The vertex geometry (NODE) has been rotated and translated
% so that the imaging plane is on the z-axis. Then we iterate
% through elements to find the containing each pixel
function EPTR= img_mapper3(NODE, ELEM, npx, npy );
  xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
  xmean= mean([xmin,xmax]); xrange= xmax-xmin;

  ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  range= max([xrange, yrange]);
  [x y]=meshgrid( ...
      linspace( xmean - range*0.5, xmean + range*0.5, npx ), ...
      linspace( ymean + range*0.5, ymean - range*0.5, npy ) );

  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A
  for j= 1: size(ELEM,2)
    xyz= NODE(:,ELEM(:,j))';
    min_z= min(xyz(:,3)); max_z= max(xyz(:,3));
    if (min_z>0 | max_z<0)
        continue;
    end
    min_x= min(xyz(:,1)); max_x= max(xyz(:,1));
    min_y= min(xyz(:,2)); max_y= max(xyz(:,2));

    % Simplex volume is det([v2-v1,v3-v1, ...])
    VOL= abs(det(xyz'*[-1,1,0,0;-1,0,1,0;-1,0,0,1]'));

    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max_y & y(:)>=min_y ...
             & x(:)<=max_x & x(:)>=min_x );

    nn=  size(ELEM,1); %Simplex vertices
    ll=  length(endr);
    vol=zeros(ll,nn);
    for i=1:nn
       i1= i; i2= rem(i,nn)+1; i3= rem(i+1,nn)+1;
       x1= xyz(i1,1)-x(endr); y1= xyz(i1,2)-y(endr); z1= xyz(i1,3);
       x2= xyz(i2,1)-x(endr); y2= xyz(i2,2)-y(endr); z2= xyz(i2,3);
       x3= xyz(i3,1)-x(endr); y3= xyz(i3,2)-y(endr); z3= xyz(i3,3);
       vol(:,i)= x1.*y2.*z3 - x1.*y3.*z2 - x2.*y1.*z3 + ...
                 x3.*y1.*z2 + x2.*y3.*z1 - x3.*y2.*z1;
    end

    endr( sum(abs(vol),2) - VOL >1e-8 )=[];
    EPTR(endr)= j;
  end %for j=1:ELEM


% Level model: usage
%   NODE= level_model( fwd_model, level );
%
% Level is a 1x3 vector specifying the x,y,z axis intercepts
% NODE describes the vertices in this coord space

function [NODE,ELEM]= level_model( fwd_model, level )
   ELEM= fwd_model.elems';

   vtx= fwd_model.nodes;
   [nn, dims] = size(vtx);
   if dims ==2 % 2D case
       NODE= vtx';
       return;
   end

   % Infinities tend to cause issues -> replace with realmax
   % Don't need to worry about the sign of the inf
   % Octave returns inf for norm of realmax
   level( isinf(level) | isnan(level) ) = 1e154; %realmax;
   level( level==0 ) =     1e-10; %eps;

   % Step 1: Choose a centre point in the plane
   %  Weight the point by it's inv axis coords
   invlev= 1./level;
   ctr= invlev / sum( invlev.^2 );

   % Step 2: Choose basis vectors in the plane
   %  First is the axis furthest from ctr
   [jnk, s_ax]= sort( - abs(level - ctr) );
   v1= [0,0,0]; v1(s_ax(1))= level(s_ax(1));
   v1= v1 - ctr;
   v1= v1 / norm(v1);

   % Step 3: Get off-plane vector, by cross product
   v2= [0,0,0]; v2(s_ax(2))= level(s_ax(2));
   v2= v2 - ctr;
   v2= v2 / norm(v2);
   v3= cross(v1,v2);

   % Step 4: Get orthonormal basis. Replace v2
   v2= cross(v1,v3);

   % Step 5: Get bases to point in 'positive directions'
   v1= v1 * (1-2*(sum(v1)<0));
   v2= v2 * (1-2*(sum(v2)<0));
   v3= v3 * (1-2*(sum(v3)<0));

   NODE= [v1;v2;v3] * (vtx' - ctr'*ones(1,nn) );
