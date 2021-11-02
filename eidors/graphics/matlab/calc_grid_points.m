function val = calc_grid_points(img, xpts, ypts, zpts)
%CALC_GRID_POINTS - image values at points in grid
% VAL = CALC_GRID_POINTS(img, xpts, ypts, zpts) returns the matrix of image
% values at points on a rectangular grid as defined by
% NDGRID(xpts,ypts,zpts).
%
% See also: GET_IMG_DATA, POINT_IN_TET, MK_GRID_MODEL, NDGRID

%TODO:
% - also output a c2f matrix
% - consider adding support for leveling the model first

% (C) 2021 Bartek Grychtol and Andy Adler. 
% License: GPL version 2 or version 3
% $Id$


if nargin == 1 && ischar(img) && strcmp(img, 'UNIT_TEST')
    do_unit_test;
    return
end

fmdl = img.fwd_model;
data = get_img_data(img); 


[fmdl, use_elem] = crop_model(fmdl, xpts, ypts, zpts);
[x,y,z]= ndgrid( xpts, ypts, zpts);
pts = [x(:),y(:),z(:)];
point2tet = point_in_tet(fmdl, pts, eps);

if isfield(img, 'elem_data')
    val =  point2tet * data(use_elem);
    val = val ./ sum(point2tet,2);
elseif isfield(img, 'node_data')
    copt.fstr = 'calc_grid_points';
    copt.cache_obj = {fmdl.nodes, fmdl.elems, pts};
    point2node = eidors_cache(@calc_point_node_interpolation, {fmdl, pts, point2tet}, copt);
    val =  point2node * data;
    pts_out = ~any(point2tet,2);
    val(pts_out) = NaN;
end


val = reshape(val, size(x));
return


%-------------------------------------------------------------------------%
% Calculate node interpolation matrix
function p2n = calc_point_node_interpolation(fmdl, pts, point2tet)
    [pts_idx,els_idx] = find(point2tet);
    
    %deal with nodes shared by multiple tets
    [pts_idx, idx] = unique(pts_idx, 'stable');
    point2tet = sparse(pts_idx, els_idx(idx), 1, size(point2tet,1), size(point2tet,2));
    
    elem_ptr = find(any(point2tet,1));
    ndims = 3;

    
    pts_map = zeros(size(pts,1),1); pts_map(pts_idx) = 1:length(pts_idx);
    pts_idx = repmat(pts_idx,1,ndims+1);

    nds_idx = zeros(size(pts_idx));
    int_ptr = zeros(size(pts_idx));

    for e= elem_ptr % loop over elements containing points
        nodes_i = fmdl.elems(e,:);
        pts_i = find(point2tet(:,e));
        nds_idx(pts_map(pts_i),:) = repmat(nodes_i,length(pts_i),1);
        
%         int_fcn = inv( [ones(1,ndims+1);fmdl.nodes(nodes_i,:)'] );
%         int_ptr(pts_map(pts_i),:) = ( int_fcn *[ones(numel(pts_i),1),pts(pts_i,:)]' )';
        
        % matlab claims this is better:
        int_ptr(pts_map(pts_i),:) = ...
            ([ones(1,ndims+1);fmdl.nodes(nodes_i,:)'] \ ...
            [ones(numel(pts_i),1),pts(pts_i,:)]')';
    end
    p2n = sparse(pts_idx,nds_idx,int_ptr, size(pts,1),size(fmdl.nodes,1));
   

%-------------------------------------------------------------------------%
% Crop parts of the model outside the grid
function [mdl, use_elem] = crop_model(mdl, varargin)
    use_elem = true([size(mdl.elems,1),1]);
    for i = 1:3
        coord = mdl.nodes(:,i);
        use_elem = use_elem & any(coord(mdl.elems) >= min(varargin{i}),2);
        use_elem = use_elem & any(coord(mdl.elems) <= max(varargin{i}),2);
    end
    mdl.elems = mdl.elems(use_elem,:);
    % is it worth removing unused nodes?
            


function do_unit_test
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl,1);
   load datacom.mat A B;
   img.elem_data(A) = 1.2;
   img.elem_data(B) = 0.8;
   xvec = -.5:.05:.5;
   yvec = -1:.05:1;
   zvec = 1:.11:2;  
   rmdl = mk_grid_model([],xvec,yvec,zvec);
   middle_point = @(v) (v(1:end-1) + v(2:end))/2;
   xpts = middle_point(xvec);
   ypts = middle_point(yvec);
   zpts = middle_point(zvec);
   
   % elem data
   val = calc_grid_points(img, xpts, ypts, zpts);
   rimg = mk_image(rmdl, val(:));
   subplot(221)
   show_fem(img);
   hold on
   show_fem(rimg);
   hold off
   subplot(222)
   show_slices(val);
   
   % node data
   img.fwd_model = fix_model(img.fwd_model,struct('node2elem',1));
   N2E = img.fwd_model.node2elem;
   nimg = rmfield(img,'elem_data');
   nimg.node_data = N2E * img.elem_data ./ sum(N2E,2);
   val = calc_grid_points(nimg, xpts, ypts, zpts);
   
   subplot(223)
   show_fem(nimg);
   hold on
   show_fem(rimg);
   hold off
   subplot(224)
   show_slices(val);
   