function map = mdl_slice_mapper( fmdl, maptype )
% MDL_SLICE_MAPPER: map pixels to FEM elements or nodes
%    map = mdl_slice_mapper( fmdl, maptype );
%
% USAGE:
% fmdl = fwd_model object
%     required fields
%   fmdl.mdl_slice_mapper.npx   - number of points in horizontal direction
%   fmdl.mdl_slice_mapper.npy   - number of points in vertical 
%    or
%   fmdl.mdl_slice_mapper.x_pts - vector of points in horizontal direction
%   fmdl.mdl_slice_mapper.y_pts - vector of points in vertical
%     x_pts starts at the left, and y_pts starts at the top, this means
%     that y_pts normally would run from max to min
%    or
%   fmdl.mdl_slice_mapper.resolution - number of points per unit
%   
%   fmdl.mdl_slice_mapper.level - any definition accepted by
%          LEVEL_MODEL_SLICE for a single slice
%   OR (deprecated)
%   fmdl.mdl_slice_mapper.centre and .rotate
%          are the centre point and rotation matrices around the point
%
% maptype
%    for 'elem' map is FEM element nearest the point
%    for 'node' map is FEM vertex nearest the point
%    for 'nodeinterp' map a npx x npy x (Nd+1) matrix such that for each point i,j
%       the nearby nodes are weighted with the corresponding element in the map(i,j).
%    for 'get_points' map contains the x an y vectors of points used for
%       for mapping as {x, y}. No mapping is performed.
%
%
% See also LEVEL_MODEL_SLICE

% (C) 2006-2022 Andy Adler and Bartek Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return; end

if nargin < 3, lev_no = 1; end

copt.log_level = 4;
params = {fmdl, maptype};
copt.cache_obj = {fmdl.nodes, fmdl.elems, fmdl.mdl_slice_mapper, maptype};
copt.fstr = 'mdl_slice_mapper';
switch maptype
   case {'elem', 'node', 'nodeinterp'}
      map = eidors_cache(@do_mdl_slice_mapper,      params, copt);
   case 'get_points'
      map = get_points(fmdl);  
   otherwise
      error('expecting maptype = elem or node');
end

function map = do_mdl_slice_mapper(fmdl, maptype)
    if isfield(fmdl.mdl_slice_mapper,'model_2d') && ...
               fmdl_model.mdl_slice_mapper.model_2d && ...
               size(fmdl.nodes,2) == 3
        fmdl.mdl_slice_mapper.nodes(:,3) = [];
    end
    [NODE, ELEM] = level_model( fmdl);
    fmdl.nodes = NODE';
    [x, y] = grid_the_space( fmdl);
   
    switch maptype
       case 'elem'
          map = mdl_elem_mapper(NODE, ELEM, x, y);
       case 'node'
          bnd = [];
          try bnd = fmdl.boundary; end
          map = mdl_node_mapper(NODE, ELEM, bnd, x, y);
       case 'nodeinterp'
          map = mdl_nodeinterp_mapper(NODE, ELEM, x, y);
    end

function elem_ptr = mdl_elem_mapper(NODE, ELEM, x, y)
   if size(NODE,1) ==2 %2D
      elem_ptr= img_mapper2( NODE, ELEM, x, y);
   else
      elem_ptr= img_mapper3( NODE, ELEM, x, y);
   end
   
function ninterp_ptr = mdl_nodeinterp_mapper(NODE, ELEM, x, y)
   ver = eidors_obj('interpreter_version');
   if ~ver.isoctave && ver.ver >= 9.004 % pointLocation was slow before
     ninterp_ptr = mdl_nodeinterp_mapper_triangulation(NODE, ELEM, x, y);
     return
   end
   elem_ptr = mdl_elem_mapper(NODE, ELEM, x, y);
   
   ndims = size(NODE,1);
   if  ndims == 2;  NODEz = []; else; NODEz= 0; end
   ninterp_ptr = zeros(length(x(:)),ndims+1); % reshape later
   
   elems = ELEM';
   for i= find( elem_ptr(:)>0 )' % look for all x,y inside elements
     nodes_i = elems(elem_ptr(i),:);
     ninterp_ptr(i,:) = ( [ones(1,ndims+1);NODE(:,nodes_i)] \ [1;x(i);y(i);NODEz] )';
   end
   ninterp_ptr = reshape( ninterp_ptr, size(x,1), size(x,2), ndims + 1);

function ninterp_ptr = mdl_nodeinterp_mapper_triangulation(NODE, ELEM, x, y)
   ndims = size(NODE,1);
   pts = [x(:),y(:)];
   if size(NODE,1) == 3, pts(:,3) = 0; end
   
   TR = triangulation(ELEM', NODE');
   [el, bc] =   TR.pointLocation(pts);
   bc(isnan(el),:) = 0;
   ninterp_ptr = reshape(bc,size(x,1), size(x,2), ndims + 1);
   
function node_ptr = mdl_node_mapper(NODE, ELEM, bnd, x, y)

   ver = eidors_obj('interpreter_version');
   if ~ver.isoctave
       node_ptr = node_mapper_triangulation( NODE, ELEM, x, y);
   else
       ndims = size(NODE,1);
       if ndims == 2
         % old code
         if isempty(bnd), bnd = find_boundary(ELEM'); end
         node_ptr= node_mapper( NODE, ELEM, bnd, x, y);
       else
         node_ptr= node_mapper_dsearchn( NODE, ELEM, x, y);
       end
   end
   
% in 3D this is somewhat faster in matlab. Perfomance in octave varies with size
function node_ptr = node_mapper_dsearchn( NODE, ELEM, x, y)
   if size(NODE,1) == 2
      node_ptr = dsearchn(NODE', ELEM', [x(:),y(:)], 0);
   else
      pts = [x(:),y(:)];
      pts(:,3) = 0;
      [NODE, ELEM, use_nodes] = limit_3dmodel_to_slice(NODE,ELEM);
      node_ptr = dsearchn(NODE', ELEM', pts , 0);
      in = node_ptr>0;
      node_ptr(in) = use_nodes(node_ptr(in));
   end
   node_ptr = reshape(node_ptr, size(x));
   
function node_ptr = node_mapper_triangulation( NODE, ELEM, x, y)
    TR = triangulation(ELEM',NODE');
    pts = [x(:), y(:)];
    if size(NODE,1) == 3
        pts(:,3) = 0;
    end
    id = TR.pointLocation(pts);
    in = ~isnan(id);
    node_ptr = zeros(size(in));
    node_ptr(in) = TR.nearestNeighbor(pts(in,:));
    node_ptr = reshape(node_ptr, size(x));
    

   
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
     NODEz2= NODE(3,:).^2;
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
        [~, ff] = min(dist2);
        NPTR(i,j) = ff;
     end
  end
  NPTR(~in)= 0; % outside

% Search through each element and find the points which
% are in that element
% EPTR is matrix npx x npy with a pointer to the
% element which contains it.
function EPTR= img_mapper2(NODE, ELEM, x, y );
  ver = eidors_obj('interpreter_version');
  if ver.isoctave
    id = tsearch(NODE(1,:),NODE(2,:), ELEM', x(:),y(:));
  else 
    TR = triangulation(ELEM',NODE');
    id = TR.pointLocation([x(:), y(:)]);
  end
  id(isnan(id)) = 0;
  EPTR = reshape(id,size(x));

% Search through each element and find the points which
% are in that element
% EPTR is matrix npx x npy with a pointer to the
% element which contains it.
function EPTR= img_mapper2_old(NODE, ELEM, x, y );
  [npy,npx] = size(x);
  v_yx= [-y(:),x(:)];
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

  ver = eidors_obj('interpreter_version');
  if ver.isoctave 
      img2d = mdl_3d_to_2d(NODE, ELEM);  
      id = tsearch(img2d.fwd_model.nodes(:,1),img2d.fwd_model.nodes(:,2), ...
                   img2d.fwd_model.elems, x(:),y(:));
    %  id = tsearchn(NODE(:,use_nodes)', map(ELEM(:,use_elem))', [x(:),y(:),zeros(numel(x),1)]);
      in = ~isnan(id);
      id(in) = img2d.elem_data(id(in));
      id(~in) = 0;
  else
      if ver.ver <  9.004 % pointLocation was slow before
          EPTR = img_mapper3_old(NODE, ELEM, x, y);
          return
      end
      TR = triangulation(ELEM', NODE');
      pts = [x(:),y(:)]; pts(:,3) = 0;
      id = pointLocation(TR, pts);
      id(isnan(id)) = 0;     
  end 
  EPTR = reshape(id,size(x));

function EPTR= img_mapper3_old(NODE, ELEM, x, y );
  [npy,npx] = size(x);

  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A
  idx = 1:size(ELEM,2);
  z = reshape(NODE(3,ELEM),size(ELEM));
  idx(min(z)>0 | max(z)<0) = [];
  for j= idx
    xyz= NODE(:,ELEM(:,j))';

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
    xendr = x(endr); yendr = y(endr);
    for i=1:nn
       i1= i; i2= rem(i,nn)+1; i3= rem(i+1,nn)+1;
       x1= xyz(i1,1)-xendr; y1= xyz(i1,2)-yendr; z1= xyz(i1,3);
       x2= xyz(i2,1)-xendr; y2= xyz(i2,2)-yendr; z2= xyz(i2,3);
       x3= xyz(i3,1)-xendr; y3= xyz(i3,2)-yendr; z3= xyz(i3,3);
       vol(:,i)= x1.*y2.*z3 - x1.*y3.*z2 - x2.*y1.*z3 + ...
                 x3.*y1.*z2 + x2.*y3.*z1 - x3.*y2.*z1;
    end

    endr( sum(abs(vol),2) - VOL >1e-8 )=[];
    EPTR(endr)= j;
  end %for j=1:ELEM
  
function [NODE, ELEM, use_nodes, use_elem] = limit_3dmodel_to_slice(NODE,ELEM)
    use_elem = 1:size(ELEM,2);
    z = reshape(NODE(3,ELEM),size(ELEM));
    rng = max(z(:)) - min(z(:));
    use_elem(min(z)>0.01*rng | max(z)<-0.01*rng) = [];
    use_nodes = unique(ELEM(:,use_elem));
    map = zeros(size(NODE,2),1); map(use_nodes) = 1:numel(use_nodes);
    NODE = NODE(:,use_nodes);
    ELEM = map(ELEM(:,use_elem));
  
% function to be used on a leveled model (slice at z=0)  
function [img2d, use_nodes] = mdl_3d_to_2d(NODE, ELEM)  
  [NODE, ELEM, use_nodes, use_elem] = limit_3dmodel_to_slice(NODE,ELEM);
  fmdl.nodes = NODE';
  fmdl.elems = ELEM';
  fmdl.type = 'fwd_model';
  fmdl.mdl_slice_mesher.interp_elems = false;
  img.fwd_model = fmdl;
  img.elem_data = use_elem;
  img.type = 'image';
  img2d = mdl_slice_mesher(img,[inf inf 0]);
  
% Level model: usage
%   NODE= level_model( fwd_model, level );
%
% Level is a 1x3 vector specifying the x,y,z axis intercepts
% NODE describes the vertices in this coord space

function [NODE, ELEM] = level_model( fwd_model )
   vtx= fwd_model.nodes;
   ELEM = fwd_model.elems';
   
   if mdl_dim(fwd_model) ==2 % 2D case
       NODE= vtx';
       return;
   end

   if     isfield(fwd_model.mdl_slice_mapper,'level')
       N_slices = level_model_slice(fwd_model.mdl_slice_mapper.level);
       if N_slices > 1
           warning(['Multiple slices defined on forward model. '...
               'Using the first slice.'])
       end
       NODE = level_model_slice(vtx, fwd_model.mdl_slice_mapper.level, 1);
   elseif isfield(fwd_model.mdl_slice_mapper,'centre')
       % just in case somebody was using that interface
       rotate = fwd_model.mdl_slice_mapper.rotate;
       centre = fwd_model.mdl_slice_mapper.centre;
       warning('EIDORS:MDL_SLICE_MAPPER:DeprecatedInterface',...
            'Specifying mdl_slice_mapper.rotate and .centre is deprecated')
       level = struct('centre', centre,'rotation_matrix',rotate);
       N_slices = level_model_slice(level);
       if N_slices > 1
           warning(['Multiple slices defined on forward model. '...
               'Using the first slice.'])
       end
       NODE = level_model_slice(vtx, level , 1);       
   else   error('mdl_slice_mapper: no field level or centre provided');
   end
   
   NODE = NODE{1}'; 
   
function pts = get_points(fwd_model);
   NODE = level_model( fwd_model );
   if isfield(fwd_model.mdl_slice_mapper,'model_2d') && ...
           fwd_model.mdl_slice_mapper.model_2d && size(NODE,1) == 3
       NODE(3,:) = [];
   end
   if size(NODE,1) ==2 %2D
      [x,y] = grid_the_space( fwd_model, 'only_get_points');
   else
      fmdl3 = fwd_model; fmdl3.nodes = NODE'; 
      [x,y] = grid_the_space( fmdl3 , 'only_get_points');
   end
   pts = {x, y};
  
   
% Create matrices x y which grid the space of NODE
function  [x,y] = grid_the_space( fmdl, flag);

  xspace = []; yspace = [];
  try 
     xspace =  fmdl.mdl_slice_mapper.x_pts;
     yspace =  fmdl.mdl_slice_mapper.y_pts;
  end

  if isempty(xspace)

     xmin = min(fmdl.nodes(:,1));    xmax = max(fmdl.nodes(:,1));
     xmean= mean([xmin,xmax]); xrange= xmax-xmin;

     ymin = min(fmdl.nodes(:,2));    ymax = max(fmdl.nodes(:,2));
     ymean= mean([ymin,ymax]); yrange= ymax-ymin;

    
     npx = []; npy = [];
     if all(isfield(fmdl.mdl_slice_mapper,{'npx','npy'}))
         npx  = fmdl.mdl_slice_mapper.npx;
         npy  = fmdl.mdl_slice_mapper.npy;
% This code broke UNIT TESTS for v3.10
%         range= max([xrange, yrange]); 
%         % for backward compatibility 
%         xspace = linspace( xmean - range*0.5, xmean + range*0.5, npx );
%         yspace = linspace( ymean + range*0.5, ymean - range*0.5, npy );
     else
         res = 1;
         try
             res =  fmdl.mdl_slice_mapper.resolution;
         end
         npx = ceil(xrange/res);
         npy = ceil(yrange/res);
     end
     if 1
       xdiv = xrange/npx; ydiv = yrange/npy;
       xspace = linspace( xmean - xrange/2 - xdiv/2, xmean + xrange*0.5 + xdiv/2, npx + 2);
       yspace = linspace( ymean + yrange/2 + ydiv/2, ymean - yrange*0.5 - ydiv/2, npy + 2);
       xspace = xspace(2:end-1);
       yspace = yspace(2:end-1);
     else
       xspace = linspace( xmean - xrange/2 , xmean + xrange*0.5 , npx);
       yspace = linspace( ymean + yrange/2 , ymean - yrange*0.5 , npy);
     end
  end
  if nargin > 1 && ischar(flag) && strcmp(flag, 'only_get_points')
      x = xspace; y = yspace;
      return
  end
  [x,y]=meshgrid( xspace, yspace );

function do_unit_test
% 2D NUMBER OF POINTS
   imdl = mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
   fmdl.nodes = 1e-15 * round(1e15*fmdl.nodes);
   fmdl.mdl_slice_mapper.level = [inf,inf,0];
   fmdl.mdl_slice_mapper.npx = 5;
   fmdl.mdl_slice_mapper.npy = 5;
   eptr = mdl_slice_mapper(fmdl,'elem');
   
   si = @(x,y) sub2ind([5,5],x, y);
   % points lie on nodes and edges, we allow any associated element
   els = {si(2,2), [25,34];
          si(2,4), [27,30];
          si(3,3), [1,2,3,4];
          si(4,2), [33,36];
          si(4,4), [28,31];
          };
   
   tst = eptr == [    0   37   38   39    0
                     46   25    5   27   42
                     47    8    1    6   41
                     48   33    7   28   40
                      0   45   44   43    0];
                      
   for i = 1:size(els,1)
    tst(els{i,1}) = ismember(eptr(els{i,1}), els{i,2});
   end
   unit_test_cmp('eptr01',tst,true(5,5));

   nptr = mdl_slice_mapper(fmdl,'node');
   test = [   0   27   28   29    0
             41    6    7    8   31
             40   13    1    9   32
             39   12   11   10   33
              0   37   36   35    0];
   unit_test_cmp('nptr01',nptr,test);

   nint = mdl_slice_mapper(fmdl,'nodeinterp');
   unit_test_cmp('nint01a',nint(2:4,2:4,1),[ 0.2627, 0.6909, 0.2627; 
                                             0.6906, 1,      0.6906;
                                             0.2627, 0.6909, 0.2627], 1e-3);

   fmdl.mdl_slice_mapper.npx = 6;
   fmdl.mdl_slice_mapper.npy = 4;
   eptr = mdl_slice_mapper(fmdl,'elem');
   res = [  0   49   38   38   52    0
           62   23    9   10   20   55
           63   24   14   13   19   54
            0   60   44   44   57    0];
   unit_test_cmp('eptr02',eptr,res);

   nptr = mdl_slice_mapper(fmdl,'node');
   res = [  0   27   15   16   29    0
           25    6    2    3    8   18
           24   12    5    4   10   19
            0   37   22   21   35    0];
   unit_test_cmp('nptr02',nptr,res);

% DIRECT POINT TESTS
   imdl = mk_common_model('a2c2',8); fmdl = imdl.fwd_model;
   fmdl.mdl_slice_mapper.level = [inf,inf,0];
   fmdl.mdl_slice_mapper.x_pts = linspace(-.95,.95,4);
   fmdl.mdl_slice_mapper.y_pts = [0,0.5];
   eptr = mdl_slice_mapper(fmdl,'elem');
   unit_test_cmp('eptr03',eptr,[ 47 8 6 41; 0 25 27 0]);

   nptr = mdl_slice_mapper(fmdl,'node');
   unit_test_cmp('nptr03',nptr,[ 40 13 9 32; 0 6 8 0]);

% 3D NPOINTS
   imdl = mk_common_model('n3r2',[16,2]); fmdl = imdl.fwd_model;
   fmdl.mdl_slice_mapper.level = [inf,inf,1]; % co-planar with faces
   fmdl.mdl_slice_mapper.npx = 4;
   fmdl.mdl_slice_mapper.npy = 4;
   eptr = mdl_slice_mapper(fmdl,'elem');
   si = @(x,y) sub2ind([4,4],x, y);
   % points lie on nodes and edges, we allow any associated element
   els = {si(1,2), [ 42,317];
          si(1,3), [ 30,305];
          si(2,1), [ 66,341];
          si(2,2), [243,518];
          si(2,3), [231,506];
          si(2,4), [  6,281];
          si(3,1), [ 78,353];
          si(3,2), [252,527];
          si(3,3), [264,539];
          si(3,4), [138,413];
          si(4,2), [102,377];
          si(4,3), [114,389];};
   tst = eptr == zeros(4); 
   for i = 1:size(els,1)
    tst(els{i,1}) = ismember(eptr(els{i,1}), els{i,2});
   end
   unit_test_cmp('eptr04',tst, true(4));
   nptr = mdl_slice_mapper(fmdl,'node');
   test = [   0   101    99     0
            103   116   113    97
            105   118   121   111
              0   107   109     0];
   unit_test_cmp('nptr04',nptr, test);

   fmdl.mdl_slice_mapper.level = [inf,0.01,inf];
   eptr = mdl_slice_mapper(fmdl,'elem'); 
   test = [621   792   777   555
           345   516   501   279
           343   515   499   277
            69   240   225     3];
   unit_test_cmp('eptr05',eptr,test);

   nptr = mdl_slice_mapper(fmdl,'node');
   test = [230   250   248   222
           167   187   185   159
           104   124   122    96
            41    61    59    33];
   unit_test_cmp('nptr05',nptr,test);

   nint = mdl_slice_mapper(fmdl,'nodeinterp');
   unit_test_cmp('nint05a',nint(2:3,2:3,4),[.1250,.1250;.8225,.8225],1e-3);
   
% Centre and Rotate
   fmdl.mdl_slice_mapper = rmfield(fmdl.mdl_slice_mapper,'level');
   fmdl.mdl_slice_mapper.centre = [0,0,0.9];
   fmdl.mdl_slice_mapper.rotate = eye(3);
   fmdl.mdl_slice_mapper.npx = 4;
   fmdl.mdl_slice_mapper.npy = 4;
   eptr = mdl_slice_mapper(fmdl,'elem');
   test = [  0    42    30     0
            66   243   229     6
            78   250   264   138
             0   102   114     0];
   unit_test_cmp('eptr06',eptr, test);

% SLOW
   imdl = mk_common_model('d3cr',[16,3]); fmdl = imdl.fwd_model;
   fmdl.nodes = 1e-15*round(1e15*fmdl.nodes);
   fmdl.mdl_slice_mapper.level = [inf,inf,1];
   fmdl.mdl_slice_mapper.npx = 64;
   fmdl.mdl_slice_mapper.npy = 64;
   t = cputime;
   eptr = mdl_slice_mapper(fmdl,'elem');
   txt = sprintf('eptr10 (t=%5.3fs)',cputime - t);
   test = [122809   122872   122872
           122873   122809   122749
           122873   122749   122749];
   unit_test_cmp(txt,eptr(10:12,10:12),test);  
   

% CHECK ORIENTATION
   imdl=mk_common_model('a2c0',16);
   img= mk_image(imdl,1); img.elem_data(26)=1.2;
   subplot(231);show_fem(img);
   subplot(232);show_slices(img);
   img.fwd_model.mdl_slice_mapper.npx= 20;
   img.fwd_model.mdl_slice_mapper.npy= 30;
   img.fwd_model.mdl_slice_mapper.level= [inf,inf,0];
   subplot(233);show_slices(img);
   img.fwd_model.mdl_slice_mapper.x_pts = [linspace(-1,1,23),.5];
   img.fwd_model.mdl_slice_mapper.y_pts = [linspace( 1,-1,34),.5];
   subplot(234);show_slices(img);
   
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl,1); vh= fwd_solve(img);
   load datacom.mat A B;
   img.elem_data(A) = 1.2;
   img.elem_data(B) = 0.8;
   img.calc_colours.transparency_thresh= 0.25;
   show_3d_slices(img);
 
   cuts = [inf, -2.5, 1.5; inf, 10, 1.5];
   subplot(235);  show_3d_slices(img ,[],[],[],cuts );
   
   cuts = [inf, inf, 0.5; 
           1e-10, 2e-10, inf;1e-10, 1e-10, inf;2e-10, 1e-10, inf;
           inf  , 1e-10, inf;1e-10, inf  , inf];
   subplot(236);  show_3d_slices(img ,[],[],[],cuts );
   