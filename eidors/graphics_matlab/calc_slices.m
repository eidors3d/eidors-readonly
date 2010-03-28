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
% PARAMETERS:
%   img.calc_slices.filter % Filter to be applied to images
%      Example:    img.calc_slices.filter = ones(3)/9
%
% rimg= np x np x I x L where np is 128 by default
% np can be adjusted by calc_colours('npoints')

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

try   np = img.calc_colours.npoints;
catch np = calc_colours('npoints');
end

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
   [elem_data, n_images] = get_img_data(img);
   rimg=zeros(np,np,n_images,num_levs);

   for lev_no = 1:num_levs
      level= levels( lev_no, 1:3 );

      rimg(:,:,:,lev_no) = calc_image_elems( elem_data, level, fwd_model, np);
   end
elseif isfield(img,'node_data')
   node_data= [img.node_data];
   if size(node_data,1)==1; node_data=node_data';end
   n_images= size(node_data,2);
   rimg=zeros(np,np,n_images,num_levs);

   for lev_no = 1:num_levs
      level= levels( lev_no, 1:3 );

      rimg(:,:,:,lev_no) = calc_image_nodes( node_data, level, fwd_model, np);
   end
else
   error('img does not have a data field');
end

% FILTER IMAGE
try   filt = img.calc_slices.filter; 
catch filt = []; end

if ~isempty(filt)
   rimg = filter_image(rimg, filt);
end

% Calculate an image by mapping it onto the node_ptr matrix
function rimg= calc_image_nodes( node_data, level, fwd_model, np);

   % elem_ptr_table also depends on the number of mapped points
   fwd_model.calc_slices.mapping_npoints=np;

   % Get node_ptr from cache, if available
   % NPtable is cell array of elem_ptrs for different levels
   NPtable = eidors_obj('get-cache', fwd_model, 'node_ptr_table');
   level_hash= eidors_var_id( level );

   if isfield(NPtable, level_hash);
      node_ptr= getfield(NPtable, level_hash);
   else
      NODE = level_model( fwd_model, level );
      node_ptr= node_mapper( NODE, fwd_model.boundary, np, np);

      NPtable = setfield(NPtable, level_hash, node_ptr);
      eidors_obj('set-cache', fwd_model, 'node_ptr_table', NPtable);
      eidors_msg('show_slices: setting cached value', 3);
   end

   backgnd= NaN;
   n_images= size(node_data,2);
   rval= [backgnd*ones(1,n_images); node_data];
   rimg= reshape( rval(node_ptr+1,:), np,np, n_images );

% Calculate an image by mapping it onto the elem_ptr matrix
function rimg= calc_image_elems( elem_data, level, fwd_model, np)

   % elem_ptr_table also depends on the number of mapped points
   fwd_model.mdl_slice_mapper.npx  = np;
   fwd_model.mdl_slice_mapper.npy  = np;
   fwd_model.mdl_slice_mapper.level= level;
   elem_ptr = mdl_slice_mapper( fwd_model, 'elem' );

   backgnd= NaN;
   n_images= size(elem_data,2);
   rval= backgnd*ones(size(elem_data)+[1,0]);
   rval(2:end,:) = elem_data;
   rimg= reshape( rval(elem_ptr+1,:), np,np, n_images );


function  rimg = filter_image(rimg, filt);
   [sz1,sz2,sz3,sz4] = size(rimg);
   for j1 = 1:sz3; for j2 = 1:sz4; 
      rsl = rimg(:,:,j1,j2);

      rna = isnan(rsl);
      rsl(rna) = 0;
      rsl = conv2(rsl, filt, 'same');
      rsl(rna) = NaN;

      rimg(:,:,j1,j2) = rsl;
   end; end

function do_unit_test
   img = calc_jacobian_bkgnd( mk_common_model('a2c2',8));
   img.calc_colours.npoints = 8; 

   imn = calc_slices(img);
   imt = NaN*ones(8); imt(3:6,2:7) = 1; imt(2:7,3:6) = 1; 
   do_indiv_test('cs 2d 1', imn, imt);

   img.elem_data(1:4) = 2;
   imn = calc_slices(img);
   imt(4:5,4:5) = 2; 
   do_indiv_test('cs 2d 2', imn, imt);

   img = calc_jacobian_bkgnd( mk_common_model('n3r2'));
   img.calc_colours.npoints = 8; 
   imn = calc_slices(img,[inf,inf,1]);

   imt = NaN*ones(8); imt(3:6,2:7) = 1; imt(2:7,3:6) = 1; 
   do_indiv_test('cs 3d 1', imn, imt);

   imn = calc_slices(img,[inf,0,inf]);
   imt = NaN*ones(8); imt(1:8,3:6) = 1; 
   do_indiv_test('cs 3d 2', imn, imt);

function do_indiv_test(txt,a,b, tol)
   if nargin < 4; tol = 0; end
   fprintf('%10s = ',txt);
   ok='fail';
   try; if isnan(a) == isnan(b); a(isnan(a))=0; b(isnan(b))=0; end; end
   try; if all(abs(a - b) <= tol);  ok='ok'; end; end
   disp(ok)
