function out_img= show_slices( img, levels )
% out_img = show_slices (img, levels ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or a array of structs
%          or output of CALC_SLICES (or equivalent multi-dimensional
%          picture array to be processed by mk_mosaic
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
%   img.show_slices.sep = number of pixels between the individual images
%   img.calc_colours.npoints = pixel width/height to map to
%
% if levels is scalar, then make levels equispaced horizontal
%          cuts through the object
%
% Show slices is roughly equivalent to:
%   rimg = calc_slices(img,levels); 
%   rimg = calc_colours(rimg,img); image(rimg);
%
% See also: CALC_SLICES, MK_MOSAIC

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

sep = 0;
try sep = img(1).show_slices.sep;
end


% TODO: because we wanted (back in 2005) to let show_slices
% handle lots of different scenarios (calling with images and
% without. It is now crufty ... and should be changed.
do_calc_slices = 0;
try if strcmp(img(1).type,'image'); do_calc_slices= 1; end;end 

if nargin<=1;
   try   levels = img(1).show_slices.levels
   catch levels = [];
   end
end

if isempty(levels) && do_calc_slices && size(img(1).fwd_model.nodes,2)==2
   levels= [Inf,Inf,0];
end

vh = [];
if size(levels,2) == 5
   vh = levels(:,4:5);
   levels = levels(:,1:3);
elseif size(levels)== [1,1]
   if size(img(1).fwd_model.nodes,2) == 2 % Can't do horiz slices for 2D model
      eidors_msg('Can''t do horizontal slices for 2D model. Showing 2D slice');
      levels= [Inf,Inf,0];
   else
      zmax= max(img(1).fwd_model.nodes(:,3));
      zmin= min(img(1).fwd_model.nodes(:,3));
      levels = linspace(zmax,zmin, levels+2);
      levels = levels(2:end-1)';
      levels = [(1+0*levels)*[Inf,Inf], levels];
   end
end

if do_calc_slices
   rimg= calc_slices( img, levels(:,1:3) );
else
   rimg= img;
end

n_col = 0;
try  n_col = img(1).show_slices.img_cols;
end

r_img = mk_mosaic(rimg, sep, vh, n_col);

c_img = calc_colours( r_img, img);
out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);


image(out_img);
axis('image');axis('off');axis('equal');

if nargout==0; clear('out_img'); end

function do_unit_test
   clf

   img=calc_jacobian_bkgnd(mk_common_model('a2t3',8)); 
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1);
   subplot(4,4,1); show_slices(img) 

   img.calc_colours.npoints= 128;
   subplot(4,4,2); show_slices(img) 

   img.calc_colours.npoints= 32;
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1:3);
   subplot(4,4,3); show_slices(img) 


   img.show_slices.img_cols= 1;
   subplot(4,4,4); show_slices(img) 

   imgn = rmfield(img,'elem_data');
   imgn.node_data=toeplitz(1:size(img.fwd_model.nodes,1),1);

   img.elem_data = img.elem_data(:,1);
   img.fwd_model.mdl_slice_mapper.npx = 10;
   img.fwd_model.mdl_slice_mapper.npy = 20;
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   subplot(4,4,5); show_slices(img);

   img.elem_data = img.elem_data(:,1);
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-100,100,20);
   img.fwd_model.mdl_slice_mapper.y_pts = linspace(-150,150,30);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   subplot(4,4,6); show_slices(img);

   subplot(4,4,7); show_slices(imgn) 

   imgn.fwd_model.mdl_slice_mapper.x_pts = linspace(-100,100,20);
   imgn.fwd_model.mdl_slice_mapper.y_pts = linspace(-150,150,30);
   imgn.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   subplot(4,4,8); show_slices(imgn) 


% 3D images
   img=calc_jacobian_bkgnd(mk_common_model('n3r2',[16,2])); 
   img.calc_colours.npoints= 16;
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1);
   subplot(4,4,9); show_slices(img,2) 


   img.elem_data=img.elem_data*[1:3];
   subplot(4,4,10); show_slices(img,2) 

   img.elem_data=img.elem_data(:,1:2);
   subplot(4,4,11); show_slices(img,[inf,inf,1;0,inf,inf;0,1,inf]);

   img.show_slices.sep = 5;
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-1,1,20);
   img.fwd_model.mdl_slice_mapper.y_pts = linspace(-1,1,30);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];

   subplot(4,4,12); show_slices(img,2) 


   img.elem_data = img.elem_data(:,1);
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,1;
           0,1,inf,3,1];
   subplot(4,4,13); show_slices(img,levels) 

   levels=[inf,inf,1,1,1;
           0,inf,inf,2,2;
           0,1,inf,  1,3];
   subplot(4,4,14); show_slices(img,levels) 

   img.elem_data = img.elem_data(:,[1,1]);
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,1;
           0,1,inf,  3,1];
   subplot(4,4,15); show_slices(img,levels) 
   
   m = calc_slices(img,levels);
   subplot(4,4,16); show_slices(m) 
   

