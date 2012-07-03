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
%   img.calc_slices.scale  % Scaling to apply to images
%   img.calc_slices.levels % an alternative way to specify levels
%
% rimg= np x np x I x L where np is 128 by default
%
%   np can be adjusted by calc_colours('npoints')
%     or by setting
%   img.fwd_model.mdl_slice_mapper.{npx, npy}
%        see help of mdl_slice_mapper for more options
%

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

np = calc_colours('npoints');
try   np = img(1).calc_colours.npoints;
end

if nargin < 2
    try 
        levels = img.calc_slices.levels; 
    catch
        levels = [];
    end
end

% Assume all fwd_models are same dimension (all 3D or 2D no mixed dims)
if mdl_dim(img(1))==2 
   if nargin>1 && ~isempty(levels);
       if ~all(levels(1,:) == [inf,inf,0])
          warning('specified levels ignored for 2D FEM');
       end
   end
   levels= [Inf,Inf,0];
elseif mdl_dim(img(1))==3 && isempty(levels)
   levels = [Inf Inf mean(img.fwd_model.nodes(:,3))];
   eidors_msg('calc_slices: no levels specified, assuming an xy plane',2);
end


rimg = [];
for i=1:length(img)
   rimg = cat(3, rimg, calc_this_slice( img(i), levels, np) );
end

function rimg = calc_this_slice( img, levels, np)
    % If scalar levels then we just create that many cut planes on z-dimension
    fwd_model = img.fwd_model;
    if size(levels)== [1,1]
       zmax= max(fwd_model.nodes(:,3));
       zmin= min(fwd_model.nodes(:,3));
       levels = linspace(zmin,zmax, levels+2);
       levels = levels(2:end-1)'*[Inf,Inf,1];
    end
    num_levs= size(levels,1);
    if isfield(img,'elem_data')
       [elem_data, n_images] = get_img_data(img);

       clear rimg;
       for lev_no = fliplr(1:num_levs)
           % start at max so memory is allocated once
          level= levels( lev_no, 1:3 );

          rimg(:,:,:,lev_no) = calc_image_elems( elem_data, level, fwd_model, np);
       end
    elseif isfield(img,'node_data')
       node_data= [img.node_data];
       if size(node_data,1)==1; node_data=node_data';end
       n_images= size(node_data,2);
       clear rimg;

       for lev_no = fliplr(1:num_levs)
          level= levels( lev_no, 1:3 );
          rimg(:,:,:,lev_no) = ...
             calc_image_nodes( node_data, level, fwd_model, np);
       end
    else
       error('img does not have a data field');
    end

    % FILTER IMAGE
    try   filt = img.calc_slices.filter; 
    catch filt = 1; end
    try   scal = img.calc_slices.scale; 
    catch scal = 1; end

    
    if filt*scal ~= 1
       filt = scal * ( filt/sum(filt(:)) );
       rimg = filter_image(rimg, filt);
    end

% Calculate an image by mapping it onto the node_ptr matrix
% This makes a blocky image to nearest node -> no longer used
function rimg= calc_image_nearestnodes( node_data, level, fwd_model, np);
   if ~isfield(fwd_model,'mdl_slice_mapper');
      fwd_model.mdl_slice_mapper.npx  = np;
      fwd_model.mdl_slice_mapper.npy  = np;
      fwd_model.mdl_slice_mapper.level= level;
   end
   node_ptr = mdl_slice_mapper( fwd_model, 'node' );

   backgnd= NaN;
   n_images= size(node_data,2);
   rval= [backgnd*ones(1,n_images); node_data];
   rimg= reshape( rval(node_ptr+1,:), [size(node_ptr), n_images]);

% Calculate an image by interpolating it onto the elem_ptr matrix
function rimg= calc_image_nodes( node_data, level, fwd_model, np)

   if ~isfield(fwd_model,'mdl_slice_mapper');
      fwd_model.mdl_slice_mapper.npx  = np;
      fwd_model.mdl_slice_mapper.npy  = np;
      fwd_model.mdl_slice_mapper.level= level;
   end

   nd_interp= mdl_slice_mapper( fwd_model, 'nodeinterp' );
   elem_ptr = mdl_slice_mapper( fwd_model, 'elem' );
   [sx,sy]= size(elem_ptr);

   node_ptr = fwd_model.elems; node_ptr = [0*node_ptr(1,:);node_ptr];
   node_ptr = reshape( node_ptr( elem_ptr+1, :), sx, sy, []);

   n_images = size(node_data,2);
   rimg= zeros(sx, sy, n_images);
   backgnd= NaN;
   
   for ni = 1:n_images
     znd = [backgnd;node_data(:,ni)]; % add NaN for background
     rimg(:,:,ni) = sum( znd(node_ptr+1) .* nd_interp, 3); 
   end


% Calculate an image by mapping it onto the elem_ptr matrix
function rimg= calc_image_elems( elem_data, level, fwd_model, np)

   if ~isfield(fwd_model,'mdl_slice_mapper');
      fwd_model.mdl_slice_mapper.npx  = np;
      fwd_model.mdl_slice_mapper.npy  = np;
      fwd_model.mdl_slice_mapper.level= level;
      % grid model sets mdl_slice_mapper.np* but not level
   elseif ~isfield(fwd_model.mdl_slice_mapper,'level');
      fwd_model.mdl_slice_mapper.level= level;
   end
   elem_ptr = mdl_slice_mapper( fwd_model, 'elem' );

   backgnd= NaN;
   n_images= size(elem_data,2);
   rval= backgnd*ones(size(elem_data)+[1,0]);
   rval(2:end,:) = elem_data;
   rimg= reshape( rval(elem_ptr+1,:), [size(elem_ptr), n_images]);


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
   img = mk_image( mk_common_model('a2c2',8));
   img.calc_colours.npoints = 8; 

   imc = calc_slices(img);
   imt = NaN*ones(8); imt(3:6,2:7) = 1; imt(2:7,3:6) = 1; 
   unit_test_cmp('cs 2d 1', imc, imt);

   imn = rmfield(img,'elem_data');
   imn.node_data = ones(size(imn.fwd_model.nodes,1),1);
   imc = calc_slices(imn);
   unit_test_cmp('cs 2d 2', imc, imt, 1e-14);

   img.elem_data(1:4) = 2;
   imc = calc_slices(img);
   imt(3:6,3:6) = 1; imt(4:5,4:5) = 2; 
   unit_test_cmp('cs 2d 3', imc, imt);

   imn.node_data(1:5) = 2;
   imc = calc_slices(imn);
   imt(3:6,3:6) = 1; imt(4:5,3:6) = 1.049020821501088;
   imt(3:6,4:5) = 1.049020821501088; imt(4:5,4:5) = 2;
   unit_test_cmp('cs 2d 4', imc, imt, 1e-14);

   imn.node_data(:) = 1; imn.node_data(1) = 4;
   imc = calc_slices(imn);
   imt(3:6,3:6) = 1; imt(4:5,4:5) = 1.575633893074693; 
   unit_test_cmp('cs 2d 5', imc, imt, 1e-14);

   imn.calc_colours.npoints = 7; 
   imc = calc_slices(imn);
   imt = NaN*ones(7); imt(2:6,2:6) = 1; imt(4,1:7)= 1; imt(1:7,4)= 1;imt(4,4) = 4; 
   unit_test_cmp('cs 2d 6', imc, imt, 1e-14);


% 3D Tests
   img = calc_jacobian_bkgnd( mk_common_model('n3r2',[16,2]));
   img.calc_colours.npoints = 8; 
   imn = calc_slices(img,[inf,inf,1]);

   imt = NaN*ones(8); imt(3:6,2:7) = 1; imt(2:7,3:6) = 1; 
   unit_test_cmp('cs 3d 1', imn, imt);

   imn = calc_slices(img,[inf,0,inf]);
   imt = NaN*ones(8); imt(1:8,3:6) = 1; 
   unit_test_cmp('cs 3d 2', imn, imt);

   % Should have no effect
   img.fwd_model.nodes(:,3) = img.fwd_model.nodes(:,3)-1;
   imn = calc_slices(img,[inf,0,inf]); 
   imt = NaN*ones(8); imt(1:8,3:6) = 1; 
   unit_test_cmp('cs 3d 3', imn, imt);


% multi image struct
   img = mk_image( mk_common_model('a2c2',8));
   img.calc_colours.npoints = 8; 
   img(2) = img;
   imc = calc_slices(img); 
   imt = NaN*ones(8); imt(3:6,2:7) = 1; imt(2:7,3:6) = 1; 
   unit_test_cmp('cs mult 1', imc, cat(3,imt,imt));

   imgb = mk_image( mk_common_model('b2c2',8));
   imgb.calc_colours.npoints = 8; 
   img(3)=imgb;
   imc = calc_slices(img); 
   unit_test_cmp('cs mult 2', imc, cat(3,imt,imt,imt));

   imdl = mk_common_model('c2t2',16);
   img = mk_image(imdl,1);
   imc = calc_slices(img);
   unit_test_cmp('size e  1', size(imc), [64,64]);

   img.calc_colours.npoints = 40;
   imc = calc_slices(img);
   unit_test_cmp('size e  2', size(imc), [40,40]);

   img.fwd_model.mdl_slice_mapper.npx = 22;
   img.fwd_model.mdl_slice_mapper.npy = 32;
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   imc = calc_slices(img);
   unit_test_cmp('size e  3', size(imc), [32,22]);

   img.fwd_model = rmfield(img.fwd_model,'mdl_slice_mapper');
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-150,150,20);
   img.fwd_model.mdl_slice_mapper.y_pts =-linspace(-150,150,23);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   imc = calc_slices(img);
   unit_test_cmp('size e  4', size(imc), [23,20]);

   img = rmfield(img,'elem_data');
   img.node_data =  ones(size(img.fwd_model.nodes,1),1);
   img.node_data =  (1:size(img.fwd_model.nodes,1))';
   img.fwd_model = rmfield(img.fwd_model,'mdl_slice_mapper');
   imc = calc_slices(img);
   unit_test_cmp('size n  1', size(imc), [40,40]);

   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-150,150,20);
   img.fwd_model.mdl_slice_mapper.y_pts =-linspace(-150,150,23);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   imc = calc_slices(img);
   unit_test_cmp('size n  2', size(imc), [23,20]);
