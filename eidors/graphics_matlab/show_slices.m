function out_img= show_slices( img, levels )
% out_img = show_slices (img, levels ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or a array of structs
% out_img= matrix in the current colormap which we can use image(out_img);
%
% PARAMETERS:
%  levels = Matrix [Lx3] of L image levels
%           each row of the matrix specifies the intercepts
%           of the slice on the x, y, z axis. To specify a z=2 plane
%           parallel to the x,y: use levels= [inf,inf,2]
%  
%  if levels is [L x 5] then levels= [x,y,z,h,v] where,
%           x,y,z specify the axes intercepts, and 
%           h,v   specify the horizontal, vertical position
%                 of that slice in the output image
% 
% IMAGE PARAMETERS:
%   img.show_slices.levels (same as above);
%   img.show_slices.img_cols = number of columns in image
%   img.calc_colours.npoints = pixel width/height to map to
%
% if levels is scalar, then make levels equispaced horizontal
%          cuts through the object
%
% Show slices is roughly equivalent to:
%   rimg = calc_slices(img,levels); 
%   rimg = calc_colours(rimg,img); image(rimg);

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

np = calc_colours('npoints');
try   np = img.calc_colours.npoints;
end

do_calc_slices = 0;
try if strcmp(img.type,'image'); do_calc_slices= 1; end;end 

if nargin<=1;
   levels= [];
end

try   levels = img.show_slices.levels
catch levels = [];
end

if isempty(levels) && do_calc_slices && size(img(1).fwd_model.nodes,2)==2
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

if do_calc_slices
   rimg= calc_slices( img, levels(:,1:3) );
else
   rimg= img;
   np = size(rimg,1);
end


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
   try   img_cols = img.show_slices.img_cols;
   catch img_cols = ceil( n_frames/vert_rows );
   end
   img_rows = ceil(n_frames/img_cols);
   img_rows = ceil(img_rows/n_levels)*n_levels; % Ensure divisible by n_levels
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
