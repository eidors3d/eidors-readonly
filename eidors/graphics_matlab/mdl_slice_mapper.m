function map = mdl_slice_mapper( fmdl, maptype );
% MDL_SLICE_MAPPER: map pixels to FEM elements or nodes
%    map = mdl_slice_mapper( fmdl, levels, maptype );
%
% USAGE:
% fmdl = fwd_model object
%     required fields
%   fmdl.mdl_slice_mapper.npx   - number of points in horizontal direction
%   fmdl.mdl_slice_mapper.npy   - number of points in vertical 
%    or
%   fmdl.mdl_slice_mapper.x_pts - vector of points in horizontal direction
%   fmdl.mdl_slice_mapper.y_pts - vector of points in vertical
%
%   fmdl.mdl_slice_mapper.level = Vector [1x3] of intercepts
%          of the slice on the x, y, z axis. To specify a z=2 plane
%          parallel to the x,y: use levels= [inf,inf,2]
% maptype
%    for 'elem' map is FEM element nearest the point
%    for 'node' map is FEM vertex nearest the point
%    for 'nodeinterp' map a npx x npy x (Nd+1) matrix such that for each point i,j
%       the nearby nodes are weighted with the corresponding element in the map(i,j).

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

switch maptype
  case 'elem';       map = mdl_elem_mapper(fmdl);
  case 'node';       map = mdl_node_mapper(fmdl);
  case 'nodeinterp'; map = mdl_nodeinterp_mapper(fmdl);
  otherwise;   error('expecting maptype = elem or node');
end

function elem_ptr = mdl_elem_mapper(fwd_model);
   level= fwd_model.mdl_slice_mapper.level;

   elem_ptr = eidors_obj('get-cache', fwd_model, 'elem_ptr');

   if ~isempty(elem_ptr)
      return;
   end

   NODE = level_model( fwd_model, level );
   ELEM= fwd_model.elems';
   if size(NODE,1) ==2 %2D
      [x,y] = grid_the_space( fwd_model);
      elem_ptr= img_mapper2( NODE, ELEM, x, y);
   else
      fmdl3 = fwd_model; fmdl3.nodes = NODE'; 
      [x,y] = grid_the_space( fmdl3 );
      elem_ptr= img_mapper3( NODE, ELEM, x, y);
   end

   eidors_obj('set-cache', fwd_model, 'elem_ptr', elem_ptr);
   eidors_msg('mdl_slice_mapper: setting cached value', 3);

function ninterp_ptr = mdl_nodeinterp_mapper(fwd_model);
   level= fwd_model.mdl_slice_mapper.level;

   ninterp_ptr = eidors_obj('get-cache', fwd_model, 'ninterp_ptr');
   if ~isempty(ninterp_ptr); return; end


   elem_ptr = mdl_elem_mapper(fwd_model);
   NODE = level_model( fwd_model, level );
   fwd_model.nodes = NODE';
   [x,y] = grid_the_space( fwd_model);

   ndims = size(NODE,1);
   if  ndims == 2;  NODEz = []; else; NODEz= 0; end
   ninterp_ptr = zeros(length(x(:)),ndims+1); % reshape later

   for i= find( elem_ptr(:)>0 )'; % look for all x,y inside elements
     nodes_i = fwd_model.elems(elem_ptr(i),:);
     int_fcn = inv( [ones(1,ndims+1);NODE(:,nodes_i)] );
     ninterp_ptr(i,:) = int_fcn *[1;x(i);y(i);NODEz];
   end
   ninterp_ptr = reshape( ninterp_ptr, size(x,1), size(x,2), ndims + 1);


   eidors_obj('set-cache', fwd_model, 'ninterp_ptr', ninterp_ptr);
   eidors_msg('mdl_slice_mapper: setting cached value', 3);

function node_ptr = mdl_node_mapper(fwd_model);
   level= fwd_model.mdl_slice_mapper.level;

   node_ptr = eidors_obj('get-cache', fwd_model, 'node_ptr');

   if ~isempty(node_ptr)
      return;
   end

   NODE = level_model( fwd_model, level );
   [x,y] = grid_the_space( fwd_model);
   node_ptr= node_mapper( NODE, fwd_model.elems', fwd_model.boundary, x, y);

   eidors_obj('set-cache', fwd_model, 'node_ptr', node_ptr);
   eidors_msg('mdl_slice_mapper: setting cached value', 3);


% Search through each element and find the points which
% are in that element
% NPTR is matrix npx x npy with a pointer to the
% node closest to it.
function NPTR= node_mapper( NODE, ELEM, bdy, x, y);
  [npy,npx] = size(x);

  NODEx= NODE(1,:);
  NODEy= NODE(2,:);
  if size(NODE,1) == 2
     NODEz2= 0;
     bdy= unique(bdy(:));
     in = inpolygon(x(:),y(:),NODE(1,bdy)',NODE(2,bdy)');
  else
     NODEz2= 0; %NODE(3,:).^2; NODEs are mapped onto the z=0 plane
     % This is a slow way to get the elems outside the space, but I don't see another
     EPTR= img_mapper3(NODE, ELEM, x, y );
     in = EPTR>0;
  end
  NPTR=zeros(npy,npx);

% This next operation can be vectorized, but we don't
%  do it because that can make really big matrices

  for i= 1: npy
     for j= 1: npx
        dist2 = (NODEx-x(i,j)).^2 + (NODEy-y(i,j)).^2 + NODEz2;
        ff = find(dist2 == min(dist2));
        NPTR(i,j) = ff(1);
     end
  end
  NPTR(~in)= 0; % outside

% Search through each element and find the points which
% are in that element
% EPTR is matrix npx x npy with a pointer to the
% element which contains it.
function EPTR= img_mapper2(NODE, ELEM, x, y );
  [npy,npx] = size(x);
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
  [x,y] = grid_the_space(npx, npy);

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
function EPTR= img_mapper3(NODE, ELEM, x, y );
  [npy,npx] = size(x);

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

function NODE= level_model( fwd_model, level )

   vtx= fwd_model.nodes;
   [nn, dims] = size(vtx);
   if dims ==2 % 2D case
       NODE= vtx';
       return;
   end

   % Infinities tend to cause issues -> replace with realmax
   % Don't need to worry about the sign of the inf
   level( isinf(level) | isnan(level) ) = realmax;
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

% Create matrices x y which grid the space of NODE
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
% 2D NUMBER OF POINTS
   imdl = mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
   fmdl.mdl_slice_mapper.level = [inf,inf,0];
   fmdl.mdl_slice_mapper.npx = 5;
   fmdl.mdl_slice_mapper.npy = 5;
   eptr = mdl_slice_mapper(fmdl,'elem');
   do_indiv_test('eptr01',eptr,[ 0  0 51  0  0; 0 34 26 30  0;
                 62 35  4 29 55; 0 36 32 31  0; 0  0 59  0  0]);

   nptr = mdl_slice_mapper(fmdl,'node');
   do_indiv_test('nptr01',nptr,[ 0  0 28  0  0; 0 14  7 17  0;
                 40 13  1  9 32; 0 23 11 20  0; 0  0 36  0  0]);

   nint = mdl_slice_mapper(fmdl,'nodeinterp');
   do_indiv_test('nint01a',nint(2:4,2:4,1),[ 0.8284, 1, 0.8284;1,1,1; 0.8284, 1, 0.8284], 1e-3);

   fmdl.mdl_slice_mapper.npx = 5;
   fmdl.mdl_slice_mapper.npy = 3;
   eptr = mdl_slice_mapper(fmdl,'elem');
   do_indiv_test('eptr02',eptr,[  0  0 51 0  0;62 35  4 29 55; 0 0 59 0 0]);

   nptr = mdl_slice_mapper(fmdl,'node');
   do_indiv_test('nptr02',nptr,[ 0 0 28 0 0; 40 13 1 9 32; 0 0 36 0 0 ]);

% DIRECT POINT TESTS
   imdl = mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
   fmdl.mdl_slice_mapper.level = [inf,inf,0];
   fmdl.mdl_slice_mapper.x_pts = linspace(-1,1,5);
   fmdl.mdl_slice_mapper.y_pts = [0,0.5];
   eptr = mdl_slice_mapper(fmdl,'elem');
   do_indiv_test('eptr03',eptr,[ 62 35 4 29 55; 0 34 26 30 0]);

   nptr = mdl_slice_mapper(fmdl,'node');
   do_indiv_test('nptr03',nptr,[ 40 13 1 9 32; 0 14 7 17 0]);

% 3D NPOINTS
   imdl = mk_common_model('n3r2',8); fmdl = imdl.fwd_model;
   fmdl.mdl_slice_mapper.level = [inf,inf,1];
   fmdl.mdl_slice_mapper.npx = 4;
   fmdl.mdl_slice_mapper.npy = 4;
   eptr = mdl_slice_mapper(fmdl,'elem');
   test = zeros(4); test(2:3,2:3) = [512 228;524 533];
   do_indiv_test('eptr04',eptr, test);
   nptr = mdl_slice_mapper(fmdl,'node');
   test = zeros(4); test(2:3,2:3) = [116 113;118 121];
   do_indiv_test('nptr04',nptr, test);

   fmdl.mdl_slice_mapper.level = [inf,0,inf];
   eptr = mdl_slice_mapper(fmdl,'elem');
   test = zeros(4); test(1:4,2:3) = [ 792 777; 791 776; 515 500; 239 224];
   do_indiv_test('eptr05',eptr,test);

   nptr = mdl_slice_mapper(fmdl,'node');
   test = zeros(4); test(1:2,:) = [ 80, 124, 122, 64; 17, 61, 59, 1];
   do_indiv_test('nptr05',nptr,test);

   nint = mdl_slice_mapper(fmdl,'nodeinterp');
   do_indiv_test('nint05a',nint(2:3,2:3,1),[0,1;0,1],1e-3);




function do_indiv_test(txt,a,b, tol)
   if nargin < 4; tol = 0; end
   fprintf('%10s = ',txt);
   ok='fail';
   try; if isnan(a) == isnan(b); a(isnan(a))=0; b(isnan(b))=0; end; end
   try; if all(abs(a - b) <= tol);  ok='ok'; end; end
   disp(ok)
