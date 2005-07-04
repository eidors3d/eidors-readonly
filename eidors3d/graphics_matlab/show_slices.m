function rimg_out = show_slices( img, levels, clim )
% show_slices (img, levels, clim  ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct
% levels = Matrix [Lx3] of L image levels
%          each row of the matrix specifies the x,y,z intercepts
%          of the slice
% clim   = colourmap limit (or default if not specified)
%        = [] => Autoscale

% $Id: show_slices.m,v 1.3 2005-07-04 09:07:09 aadler Exp $

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

function rimg= calc_image( img, level, clim)

fwd_model= img.fwd_model;
elem_ptr = eidors_obj('get-cache', fwd_model, 'elem_ptr');

np= 128;

if ~isempty(elem_ptr)
   eidors_msg('show_slices: using cached value', 2);
else
   NODE= fwd_model.nodes';
   ELEM= fwd_model.elems';
   elem_ptr= img_mapper( NODE, ELEM, np, np);
   eidors_obj('set-cache', fwd_model, 'elem_ptr', elem_ptr);
end




backgnd= .01;
scale=1;
rval= [backgnd; scale*img.elem_data];
%rimg= reshape( rval(eptr+1), npy,npx );
rimg= reshape( rval(elem_ptr+1), np,np );


% create a set of pointers to elements in a 2D representation
% of the mesh
function EPTR= img_mapper(NODE, ELEM, npx, npy );
  [x y]=meshgrid( ...
      linspace( min(NODE(1,:))*1.05, max(NODE(1,:))*1.05 ,npx ), ...
     -linspace( min(NODE(2,:))*1.05, max(NODE(2,:))*1.05 ,npy )  ); 
  v_yx= [-y(:) x(:)];
  tourne= [0 -1 1;1 0 -1;-1 1 0];
  EPTR=zeros(npy,npx);
  for j= 1: size(ELEM,2)
    xy= NODE(:,ELEM(:,j))';
    a= xy([2;3;1],1).*xy([3;1;2],2)- xy([3;1;2],1).*xy([2;3;1],2);
    endr=find( y(:)<=max(xy(:,2)) & y(:)>=min(xy(:,2)) ...
             & x(:)<=max(xy(:,1)) & x(:)>=min(xy(:,1)) );
    aa= sum(abs(ones(length(endr),1)*a'+ ...
                v_yx(endr,:)*xy'*tourne)');
    endr( abs( (abs(sum(a))-aa) ./ sum(a)) >1e-8)=[];
    EPTR(endr)= j;
  end %for j=1:ELEM

