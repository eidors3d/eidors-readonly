function out_img= show_slices( img, levels )
% out_img = show_slices (img, levels ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or a array of structs
% out_img= matrix in the current colormap which we can use image(out_img);
%
% levels = Matrix [Lx3] of L image levels
%          each row of the matrix specifies the intercepts
%          of the slice on the x, y, z axis. To specify a z=2 plane
%          parallel to the x,y: use levels= [inf,inf,2]
% 
% if levels is [L x 5] then levels= [x,y,z,h,v] where,
%          x,y,z specify the axes intercepts, and 
%          h,v   specify the horizontal, vertical position
%                of that slice in the output image
%
% if levels is scalar, then make levels equispaced horizontal
%          cuts through the object

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: show_slices.m,v 1.51 2008-07-19 17:45:08 aadler Exp $

try   np = img.calc_colours.npoints;
else  np = calc_colours('npoints');
end

dims= size(img(1).fwd_model.nodes,2);
if nargin<=1;
   levels= [];
end

if isempty(levels) && dims==2
   levels= [Inf,Inf,0];
end

if size(levels,2) == 5
   spec_position= 1;
elseif size(levels)== [1,1]
   zmax= max(img(1).fwd_model.nodes(:,3));
   zmin= min(img(1).fwd_model.nodes(:,3));
   levels = linspace(zmax,zmin, levels+2);
   levels = levels(2:end-1)';
   levels = [(1+0*levels)*[Inf,Inf], levels];
   spec_position= 0;
else
   spec_position= 0;
end

rimg= calc_slices( img, levels(:,1:3) );

n_frames = size(rimg,3);
n_levels = size(rimg,4);
if spec_position %won't work for multiple image inputs
   img_cols = max( levels(:,4) );
   img_rows = max( levels(:,5) );
else
   % vertical slices must be kept together
   % To nicely fill the image: img_cols ~ img_rows
   % Thus,  n_frames/vert_rows ~ vert_rows*n_levels;
   % or     vert_rows^2 ~ n_frames / n_levels
   vert_rows = ceil( sqrt(n_frames / n_levels) );
   
   img_cols = ceil( n_frames/vert_rows );
   img_rows = vert_rows*n_levels;
end

r_img = NaN*ones(img_rows*np, img_cols*np);

idx= (-np:-1)+1;
imno= 1;
for img_no = 1:n_frames
   for lev_no = 1:n_levels
      if spec_position %won't work for multiple image inputs
         i_col= levels( lev_no, 4) + img_no -1;
         i_row= levels( lev_no, 5);
      else
         i_col= rem( img_no-1, img_cols) + 1;
         i_row= (ceil( img_no / img_cols) -1) * n_levels + lev_no ;
      end
% disp([imno, vert_rows, img_cols, img_rows, img_no, lev_no, i_col, i_row]);
      r_img(i_row*np + idx, i_col*np + idx) = rimg(:,:,img_no,lev_no);
      imno= imno + 1;
   end
end

c_img = calc_colours( r_img, img);
out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);


image(out_img);
axis('image');axis('off');axis('equal');

if nargout==0; clear('out_img'); end
