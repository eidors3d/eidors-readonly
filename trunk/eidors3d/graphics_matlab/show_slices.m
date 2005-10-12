function rimg_out = show_slices( img, levels, clim )
% show_slices (img, levels, clim  ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct
% levels = Matrix [Lx3] of L image levels
%          each row of the matrix specifies the intercepts
%          of the slice on the x, y, z axis. To specify a z=2 plane
%          parallel to the x,y: use levels= [inf,inf,2]
% clim   = colourmap limit (or default if not specified)
%        = [] => Autoscale

% $Id: show_slices.m,v 1.8 2005-10-12 03:12:54 aadler Exp $

% NOTES:
%  - currently works for 2D samples only
%  - 

dims= size(img.fwd_model.nodes,2);
if ~exist('levels') && dims==2
   levels= [Inf,Inf,0];
end

if ~exist('clim')
   clim = [];
end

for img_no = 1:prod(size( img ))
   for lev_no = 1:size( levels,1 )
      level= levels( lev_no, : );
      rimg= calc_image( img( img_no ), level, clim );
   end
end

if nargout==0
   imagesc( rimg );
else
   rimg_out = rimg;
end

% Calculate an image by mapping it onto the elem_ptr matrix
function rimg= calc_image( img, level, clim)

fwd_model= img.fwd_model;
elem_ptr = eidors_obj('get-cache', fwd_model, 'elem_ptr');

np= 128;

if ~isempty(elem_ptr)
   eidors_msg('show_slices: using cached value', 2);
else
   NODE= fwd_model.nodes';
   ELEM= fwd_model.elems';
   elem_ptr= img_mapper2a( NODE, ELEM, np, np);
   eidors_obj('set-cache', fwd_model, 'elem_ptr', elem_ptr);
end




backgnd= .01;
scale=1;
rval= [backgnd; scale*img.elem_data];
%rimg= reshape( rval(eptr+1), npy,npx );
rimg= reshape( rval(elem_ptr+1), np,np );


% Search through each element and find the points which
% are in that element
function EPTR= img_mapper2(NODE, ELEM, npx, npy );
  xmin = min(NODE(1,:));    xmax = max(NODE(1,:));
  xmean= mean([xmin,xmax]); xrange= xmax-xmin;

  ymin = min(NODE(2,:));    ymax = max(NODE(2,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  [x y]=meshgrid( ...
      linspace( xmean- xrange*0.55, xmean + xrange*0.55, npx ), ...
      linspace( ymean- yrange*0.55, ymean + yrange*0.55, npy ) );
  v_yx= [-y(:) x(:)];
  turn= [0 -1 1;1 0 -1;-1 1 0];
  EPTR=zeros(npy,npx);
  % for each element j, we get points on the simplex a,b,c
  %   area A = abc
  %   for each candidate point d,
  %      area AA = abd + acd + bcd
  %      d is in j if AA = A  
  for j= 1: size(ELEM,2)
    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
             & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) );
    % calculate area of three subtrianges to each candidate point.
    xy= NODE(:,ELEM(:,j))';
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

  ymin = min(NODE(1,:));    ymax = max(NODE(1,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  [x y]=meshgrid( ...
      linspace( xmean - xrange*0.55, xmean + xrange*0.55, npx ), ...
      linspace( ymean - yrange*0.55, ymean + yrange*0.55, npy ) );

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
    VOL= det(xyz'*[-1,1,0;-1,0,1]');

    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max_y & y(:)>=min_y ...
             & x(:)<=max_x & x(:)>=min_x );

    nn=  size(ELEM,1); %Simplex vertices
    ll=  length(endr);
    ot= ones(1,ll);
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

  ymin = min(NODE(1,:));    ymax = max(NODE(1,:));
  ymean= mean([ymin,ymax]); yrange= ymax-ymin;

  [x y]=meshgrid( ...
      linspace( xmean- xrange*0.55, xmean + xrange*0.55, npx ), ...
      linspace( ymean- yrange*0.55, ymean + yrange*0.55, npy ) );

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

    % The simplex volume is 1/N! of the enclosing rectange
    % However, we can skip the 1/6 (3D) factor since it is common
    VOL= (max_x-min_x)*(max_y-min_y)*(max_z-min_z);

    % come up with a limited set of candidate points which
    % may be within the simplex
    endr=find( y(:)<=max_y & y(:)>=min_y ...
             & x(:)<=max_x & x(:)>=min_x );

    ot= ones(size(endr));
    vol=zeros(size(endr));
    for i=1: size(ELEM,1) % Simplex vertices
       vv=1:4; vv(i)=[]; %eliminate this vertex
       xx= [x(endr)';xyz(vv,1)*ot]; 
       yy= [y(endr)';xyz(vv,2)*ot];
       zz= [z(endr)';xyz(vv,3)*ot];
       vol(i,:)= (max(xx)-min(xx)).* ...
                 (max(yy)-min(yy)).* ...
                 (max(zz)-min(zz));
    end
       
    endr( sum(vol) - VOL >1e-8 )=[];
    EPTR(endr)= j;
  end %for j=1:ELEM

