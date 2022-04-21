function rimg = calc_slices( img, levels );
% calc_slices (img, levels, clim  ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or array of I images
% levels = any level definition accepted by level_model_slice
%
% PARAMETERS:
%   img.calc_slices.filter % Filter to be applied to images
%      Example:    img.calc_slices.filter = ones(3)/9
%   img.calc_slices.scale  % Scaling to apply to images
%   img.get_img_data.frame_select = which frames of image to display
%
% rimg= np x np x I x L where np is 128 by default
%
%   np can be adjusted by calc_colours('npoints')
%     or by setting
%   img.fwd_model.mdl_slice_mapper.{npx, npy}
%        see help of mdl_slice_mapper for more options
%
% See also: LEVEL_MODEL_SLICE

% (C) 2006-2022 Andy Adler and Bartek Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

np = calc_colours('npoints');
try np = img(1).calc_colours.npoints; end

if nargin < 2, levels = []; end

if isfield(img, 'calc_slices') && isfield(img.calc_slices, 'levels')
    warning('EIDORS:CALC_SLICES:DeprecatedInterface', ...
        ['Ingoring deprecated img.calc_slices.levels definition. '...
        'Use img.fwd_model.mdl_slice_mapper.level instead. '...
        'See help level_model_slice for valid inputs.'])
    img.calc_slices = rmfield(img.calc_slices, 'levels');
end



img = data_mapper(img);

% Assume all fwd_models are same dimension (all 3D or 2D no mixed dims)
if mdl_dim(img(1))==2 
    if nargin>1
        eidors_msg('Specified levels ignored for 2D FEM', 4);
    end
    if isfield(img(1).fwd_model,'mdl_slice_mapper') && isfield(img(1).fwd_model.mdl_slice_mapper, 'level')
        eidors_msg('mdl_slice_mapper.level definition ignored for 2D FEM',4); 
        for i = 1:numel(img)
            img(i).fwd_model.mdl_slice_mapper = rmfield(img(i).fwd_model.mdl_slice_mapper,'level');
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
    
    if ~isfield(fwd_model,'mdl_slice_mapper')
        fwd_model.mdl_slice_mapper.npoints = np;
        fwd_model.mdl_slice_mapper.level= levels;
        % grid model sets mdl_slice_mapper.np* but not level
    elseif ~isfield(fwd_model.mdl_slice_mapper,'level')
        fwd_model.mdl_slice_mapper.level= levels;
    end
    
    nodes = {fwd_model.nodes};
    if size(fwd_model.nodes,2)==3 && size(fwd_model.elems,2) > 3 
        nodes = level_model_slice(fwd_model, levels);
        
        % we'll level by replacing nodes. Disable.
        lvl = struct('rotation_matrix',eye(3), 'centre', zeros(1,3));
        fwd_model.mdl_slice_mapper.level = lvl;
    end
    
    [data, n_images] = get_img_data(img);
    
    if ~any(isfield(img, {'elem_data', 'node_data'}))
       error('img does not have a data field');
    end
    n_levels = numel(nodes);
    % levels are numbered from min to max, but we want the max to end up on
    % top when shown via show_slices
    for lev_no = 1:n_levels
        fwd_model.nodes = nodes{lev_no};
        if isfield(img,'elem_data')
          rimg(:,:,:,n_levels + 1 - lev_no) = calc_image_elems( data, fwd_model);
        elseif isfield(img,'node_data')
          rimg(:,:,:,n_levels + 1 - lev_no) = calc_image_nodes( data, fwd_model);
        end
    end
    
    % FILTER IMAGE
    try   
        filt = img.calc_slices.filter; 
    catch
        filt = 1; 
    end
    try   
        scal = img.calc_slices.scale; 
    catch
        scal = 1; 
    end

    filt = scal * ( filt/sum(filt(:)) );
    rimg = filter_image(rimg, filt);

% Calculate an image by mapping it onto the node_ptr matrix
% This makes a blocky image to nearest node -> no longer used
function rimg= calc_image_nearestnodes( node_data, fwd_model)

   node_ptr = mdl_slice_mapper( fwd_model, 'node' );

   backgnd= NaN;
   n_images= size(node_data,2);
   rval= [backgnd*ones(1,n_images); node_data];
   rimg= reshape( rval(node_ptr+1,:), [size(node_ptr), n_images]);

% Calculate an image by interpolating it onto the elem_ptr matrix
function rimg= calc_image_nodes( node_data, fwd_model)

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
function rimg= calc_image_elems( elem_data, fwd_model)

   elem_ptr = mdl_slice_mapper( fwd_model, 'elem');

   backgnd= NaN;
   n_images= size(elem_data,2);
   rval= backgnd*ones(size(elem_data)+[1,0]);
   rval(2:end,:) = elem_data;
   rimg= reshape( rval(elem_ptr+1,:), [size(elem_ptr), n_images]);


function  rimg = filter_image(rimg, filt);
    
   %%% Total MATLAB BS )(*&#$)(*#&@
   %%% the && operator used to short circuit. Now it doesn't
   %%% How the (*&)(*& can anyone take this language seriously
   % all(size(filt*scal) == [1,1]) && filt*scal == 1
   if all(size(filt)==1) if filt == 1; return; end ; end

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
   imt = NaN*ones(8); imt(2:7,2:7) = 1; imt(3:6,:) = 1; imt(:,3:6) = 1; 
   unit_test_cmp('cs 2d 1', imc, imt);

   imn = rmfield(img,'elem_data');
   imn.node_data = ones(size(imn.fwd_model.nodes,1),1);
   imc = calc_slices(imn);
   unit_test_cmp('cs 2d 2', imc, imt, 1e-14);

   img.elem_data(1:4) = 2;
   imc = calc_slices(img);
   imt(4:5,4:5) = 2; 
   unit_test_cmp('cs 2d 3', imc, imt);

   imn.node_data(1:5) = 2;
   imc = calc_slices(imn);
   imt(3:6,3:6) = 1; imt(4:5,3:6) = 1.292893218813452;
   imt(3:6,4:5) = 1.292893218813452; imt(4:5,4:5) = 2;
   unit_test_cmp('cs 2d 4', imc, imt, 1e-14);

   imn.node_data(:) = 1; imn.node_data(1) = 4;
   imc = calc_slices(imn);
   imt(3:6,3:6) = 1; imt(4:5,4:5) = 1.878679656440358; 
   unit_test_cmp('cs 2d 5', imc, imt, 1e-14);

   imn.calc_colours.npoints = 7; 
   imc = calc_slices(imn);
   imt = NaN*ones(7); imt(3:5,:) = 1; imt(:,3:5)= 1; imt(2:6,2:6)= 1;imt(4,4) = 4; 
   unit_test_cmp('cs 2d 6', imc, imt, 1e-14);


% 3D Tests
   img = calc_jacobian_bkgnd( mk_common_model('n3r2',[16,2]));
   img.calc_colours.npoints = 8; 
   imn = calc_slices(img,[inf,inf,1]);

   imt = NaN*ones(8); imt(3:6,:) = 1; imt(:,3:6) = 1; imt(2:7,2:7) = 1; 
   unit_test_cmp('cs 3d 1', imn, imt);

   img.calc_colours.npoints = 12; 
   imn = calc_slices(img,[inf,0,inf]);
   imt = NaN*ones(12); imt(:,3:10) = 1; 
   unit_test_cmp('cs 3d 2', imn, imt);

   % Should have no effect
   img.fwd_model.nodes(:,3) = img.fwd_model.nodes(:,3)-1;
   imn = calc_slices(img,[inf,0,inf]); 
   unit_test_cmp('cs 3d 3', imn, imt);


% multi image struct
   img = mk_image( mk_common_model('a2c2',8));
   img.calc_colours.npoints = 8; 
   img(2) = img;
   imc = calc_slices(img); 
   imt = NaN*ones(8); imt(3:6,:) = 1; imt(:,3:6) = 1; imt(2:7,2:7) = 1; 
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
   imc = calc_slices(img);
   unit_test_cmp('size e  3', size(imc), [32,22]);

   img.fwd_model = rmfield(img.fwd_model,'mdl_slice_mapper');
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-150,150,20);
   img.fwd_model.mdl_slice_mapper.y_pts =-linspace(-150,150,23);
   imc = calc_slices(img);
   unit_test_cmp('size e  4', size(imc), [23,20]);

   img = rmfield(img,'elem_data');
   img.node_data =  ones(size(img.fwd_model.nodes,1),1);
   img.node_data =  (1:size(img.fwd_model.nodes,1))';
   img.fwd_model = rmfield(img.fwd_model,'mdl_slice_mapper');
   imc = calc_slices(img);
   unit_test_cmp('size n  1', size(imc), [40,40]);

   im2= img;
   im2.node_data =  ones(size(img.fwd_model.nodes,1),5);
   imc = calc_slices(im2);
   unit_test_cmp('size n x5', size(imc), [40,40,5]);
   im2.get_img_data.frame_select = 2:3;
   imc = calc_slices(im2);
   unit_test_cmp('size n x2', size(imc), [40,40,2]);

   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-150,150,20);
   img.fwd_model.mdl_slice_mapper.y_pts =-linspace(-150,150,23);
   imc = calc_slices(img);
   unit_test_cmp('size n  2', size(imc), [23,20]);
