function [nimg out] = mdl_slice_mesher(fmdl,level,varargin)
%MDL_SLICE_MESHER A slice of a 3D FEM as a 2D FEM 
% img2d = mdl_slice_mesher(mdl3d,level) returns a 2D FEM model MDL2D 
% suitable for viewing with SHOW_FEM representing a cut through MDL3D at
% LEVEL. 
% Note that where the intersection of an element of MDL3D and the LEVEL is
% a quadrangle, this will be represented as two triangles in IMG2D. 
% Faces of MDL3D co-planar with LEVEL will be assigned an avarage of the
% values of the two elements that share them. 
%
% [img2d ptch] = mdl_slice_mesher(...) also returns a struct PTCH suitable
% for use with the PATCH function. It offers the advantage of displaying
% both quad and tri elements. Colors can be controlled by adding
% additional arguments passed to calc_colours:
% [img2d ptch] = mdl_slice_mesher(mdl3d,level, ... )
% 
% Inputs:
%   MDL3D  - an EIDORS fwd_model or img struct with elem_data
%   LEVEL  - Vector [1x3] of intercepts
%          of the slice on the x, y, z axis. To specify a z=2 plane
%          parallel to the x,y: use levels= [inf,inf,2]
%
% To control the transparency use transparency_tresh (see CALC_COLOURS for
% details), e.g.:
%    img2d.calc_colours.transparency_thresh = -1; (no transperency)
%    calc_colours('transparency_thresh', 0.25); (some transparency)
%
% See also: SHOW_FEM, MDL_SLICE_MAPPER, SHOW_3D_SLICES, CROP_MODEL,
%           CALC_COLOURS. PATCH

% (C) 2012 Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

% TODO: 
%  1. More intuitive cut plane specification
%  2. Support node_data


if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return, end;

if isempty(varargin)
   try
      varargin{1}.calc_colours = fmdl.calc_colours;
   end
end

switch fmdl.type
    case 'image'  
       img = fmdl;
       fmdl = fmdl.fwd_model;
    case 'fwd_model'
       img = mk_image(fmdl,1);
    otherwise; error('Unknown object type');
end

opt.cache_obj = {fmdl.nodes, fmdl.elems, fmdl.electrode, level};
opt.fstr      = 'mdl_slice_mesher';


[nmdl f2c p_struct] = eidors_cache(@do_mdl_slice_mesher,{fmdl, level},opt);
nimg = build_image(nmdl, f2c, img);

switch nargout
   case 2
      out = draw_patch(p_struct, nimg.fwd_model, img.elem_data, varargin{:});
   case 0
      out = draw_patch(p_struct, nimg.fwd_model, img.elem_data, varargin{:});
      cmap_type = calc_colours('cmap_type');
      try 
         calc_colours('cmap_type',varargin{1}.calc_colours.cmap_type);
      end
      colormap(calc_colours('colourmap'));
      patch(out);
      calc_colours('cmap_type',cmap_type);
      clear nimg;
end


function [nmdl f2c out] = do_mdl_slice_mesher(fmdl,level)

mdl = fmdl;
opt.edge2elem = true;
opt.node2elem = true;
mdl = fix_model(mdl,opt);
edges = mdl.edges;
edge2elem = mdl.edge2elem;
tmp = mdl;
tmp.nodes = level_model( tmp, level )';
[nodeval nodedist] = nodes_above_or_below(tmp,0);
% find which edges are on electrodes
e_nodes = zeros(length(mdl.nodes),1);
try
   for i = 1:length(mdl.electrode)
      e_nodes(mdl.electrode(i).nodes) = i;
   end
end
e_edges = (sum(e_nodes(edges),2)/ 2)  .* (e_nodes(edges(:,1)) == e_nodes(edges(:,2)));
%% crossed edges
% exclude edges on plane (dealt with later)
idx = (sum(nodeval(edges),2) == 0) & (nodeval(edges(:,1)) ~= 0) ; 
dist = (nodedist(edges(idx,2)) - nodedist(edges(idx,1)));
t = -nodedist(edges(idx,1))./dist;
% new nodes along the edges
nodes = mdl.nodes(edges(idx,1),:) + ...
    repmat(t,1,3).*(mdl.nodes(edges(idx,2),:) - mdl.nodes(edges(idx,1),:));
% nn indexes the just-created nodes, els indexes elements
if any(idx)
    [nn els] = find(edge2elem(idx,:));
else
    nn = []; els = [];
end
els_edge = els;

electrode_node = e_edges(idx);
%% crossed nodes
idx = nodeval == 0;
ln = length(nodes); %store the size
nodes = [nodes; mdl.nodes(idx,:)]; % add the crossed nodes to the new model
electrode_node = [electrode_node; e_nodes(idx)];
% nnn indexes mdl.nodes(idx), eee indexes mdl.elems
[nnn eee] = find(mdl.node2elem(idx,:));
nnn = nnn + ln; %make a proper index into nodes
% n1: eee = ueee(n1)
[ueee jnk n1] = unique(eee,'last');
nodes_per_elem = jnk;
nodes_per_elem(2:end) = diff(jnk);
% if an elem has 3 crossed nodes, there must be 2 of them, add both for now
add = find(nodes_per_elem == 3);
if ~isempty(add)
    % what to do with faces shared between elements?
    addel = ueee(add*[1 1 1])';
    els = [els; addel(:)];
    for i = 1:length(add)
       addnd = nnn(n1 == add(i));
       nn = [nn; addnd(:)];
    end
    [els idx] = sort(els);
    nn = nn(idx);
end
% for elems with less than 4 crossed edges -> add crossed nodes if needed
[uels jnk n] = unique(els_edge,'last');

% only consider elements who have both a crossed node and edge
[idx ia ib] = intersect(ueee, uels);
for i = 1:length(ia)
    newnodes = nnn(n1==ia(i));
    nn = [nn; newnodes];
    els = [els; repmat(uels(ib(i)),length(newnodes),1)];
end
[els idx] = sort(els);
nn = nn(idx);
[uels jnk n] = unique(els,'last');
nodes_per_elem = jnk;
nodes_per_elem(2:end) = diff(jnk);

n_tri = length(uels) + sum(nodes_per_elem==4);
nmdl.type = 'fwd_model';
nmdl.nodes = nodes;
nmdl.elems = zeros(n_tri,3);
n_el_data = size(fmdl.elems,1);
f2c = sparse(n_el_data,length(uels));
c = 1;
% TODO: Speed this up
for i = 1:length(uels)
    switch nodes_per_elem(i)
        case 3
            nmdl.elems(c,:) = nn(n==i);
            f2c(uels(i),c) = 1;
            c = c + 1;
        case 4
            nds = nn(n==i);
            nmdl.elems(c,:) = nds(1:3);
            f2c(uels(i),c) = 1;
            nmdl.elems(c+1,:) = nds(2:4);
            f2c(uels(i),c+1) = 1;
            c = c + 2;
    end
end
% deal with double elements (from shared faces)
nmdl.elems = sort(nmdl.elems,2);
[nmdl.elems idx] = sortrows(nmdl.elems);
f2c = f2c(:,idx);
[nmdl.elems n idx] = unique(nmdl.elems, 'rows');
[x y] = find(f2c);
% put all source elements that contribute to destination element on one
% column
f2c = sparse(x,idx,1);
% ensure correct number of columns (happens when the last source element
% doesn't contribute 
if size(f2c,1) < n_el_data
   f2c(n_el_data,end) = 0;
end
n_src_els = sum(f2c,1);
f2c = f2c * spdiag(1./n_src_els);


% add electrodes
try
   for i = 1:length(mdl.electrode)
      nmdl.electrode(i) = mdl.electrode(i);
      nmdl.electrode(i).nodes = find(electrode_node == i);
   end
end

out.uels = uels;
out.els  = els;
out.nn   = nn;


function nimg = build_image(nmdl, f2c, img)
nimg = mk_image(nmdl,1);
nimg.elem_data = (img.elem_data' * f2c)';
try
   nimg.calc_colours = img.calc_colours;
end


function out = draw_patch(in, nmdl, elem_data, varargin)

uels = in.uels;
els  = in.els;
nn   = in.nn;
nodes = nmdl.nodes;

img.elem_data = elem_data;

out.Vertices = nodes;
out.Faces    = NaN(length(uels),4);
ed = zeros(length(uels),1);

% show_fem(nimg);
for i = 1:length(uels)
    idx = els == uels(i);
    id = find(idx);
    nnn = nodes(nn(idx),:);
    if nnz(idx)==4
        if intersection_test(nnn(1,:),nnn(4,:),nnn(2,:),nnn(3,:))
            id(3:4) = id([4 3]);
        elseif intersection_test(nnn(1,:),nnn(2,:),nnn(3,:),nnn(4,:))
            id(2:3) = id([3 2]);
        end
    end
%     patch(nodes(nn(id),1),nodes(nn(id),2),nodes(nn(id),3),img.elem_data(uels(i)));
    out.Faces(i,1:length(id)) = nn(id);
    ed(i) = img.elem_data(uels(i));
end
[out.FaceVertexCData scl_data] = calc_colours(ed,varargin{:});
try
   out.FaceVertexAlphaData = double(abs(scl_data) > varargin{1}.calc_colours.transparency_thresh);
   out.FaceAlpha = 'flat';
end
out.FaceColor = 'flat';
out.CDataMapping = 'direct';
% colormap(calc_colours('colourmap'));


function res = intersection_test(A,B,C,D)
% checks for intersection of segments AB and CD
% if AB and CD intersect in 3D, then their projections on a 2D plane also
% intersect (or are colliniar).

% assume they don't
res = false;

% check for interesection on the 3 cartesian planes
idx = 1:3;
for i = 0:2
    id = circshift(idx',i)';
    id = id(1:2);
    res = res || ( sign(signed_area(A(id), B(id), C(id))) ~= ...
                   sign(signed_area(A(id), B(id), D(id)))        );
end

function a = signed_area(A,B,C)
    a = ( B(1) - A(1) ) * ( C(2) - A(2) ) - ...
        ( C(1) - A(1) ) * ( B(2) - A(2) );


function [nodeval dist] = nodes_above_or_below(mdl,level)

% Set a model-dependent tolerance
tol = min(max(mdl.nodes) - min(mdl.nodes)) * 1e-10;

dist = mdl.nodes(:,3) - level;
dist(abs(dist) < tol) = 0;
nodeval = sign(dist);




% Level model: usage
%   NODE= level_model( fwd_model, level );
%
% Level is a 1x3 vector specifying the x,y,z axis intercepts
% NODE describes the vertices in this coord space

function NODE= level_model( fwd_model, level )

   vtx= fwd_model.nodes;
   [nn, dims] = size(vtx);
   if dims ==2 % 2D case
       NODE= vtx';
       return;
   end

   % Infinities tend to cause issues -> replace with realmax
   % Don't need to worry about the sign of the inf
   level( isinf(level) | isnan(level) ) = realmax;
   level( level==0 ) =     1e-10; %eps;

   % Step 1: Choose a centre point in the plane
   %  Weight the point by it's inv axis coords
   invlev= 1./level;
   ctr= invlev / sum( invlev.^2 );

   % Step 2: Choose basis vectors in the plane
   %  First is the axis furthest from ctr
   [jnk, s_ax]= sort( - abs(level - ctr) );
   v1= [0,0,0]; v1(s_ax(1))= level(s_ax(1));
   v1= v1 - ctr;
   v1= v1 / norm(v1);

   % Step 3: Get off-plane vector, by cross product
   v2= [0,0,0]; v2(s_ax(2))= level(s_ax(2));
   v2= v2 - ctr;
   v2= v2 / norm(v2);
   v3= cross(v1,v2);

   % Step 4: Get orthonormal basis. Replace v2
   v2= cross(v1,v3);

   % Step 5: Get bases to point in 'positive directions'
   v1= v1 * (1-2*(sum(v1)<0));
   v2= v2 * (1-2*(sum(v2)<0));
   v3= v3 * (1-2*(sum(v3)<0));

   NODE= [v1;v2;v3] * (vtx' - ctr'*ones(1,nn) );

   
function do_unit_test
    imdl = mk_common_model('n3r2',[16,2]);
    img = mk_image(imdl.fwd_model,1);
    load datacom.mat A B;
    img.elem_data(A) = 1.2;
    img.elem_data(B) = 0.8;
    subplot(131)
    show_fem(img);
    subplot(132)
    cla
%     show_fem(img.fwd_model);
    hold on
%     slc = mdl_slice_mesher(img, [3 -3 2]);
%     slc.calc_colours.transparency_thresh = -1;
%     show_fem(slc);
    slc = mdl_slice_mesher(img, [0 inf inf]);
    slc.calc_colours.transparency_thresh = -1;
    slc.fwd_model.boundary = slc.fwd_model.elems;
    show_fem(slc);
    slc = mdl_slice_mesher(img, [inf inf 2]);
    slc.calc_colours.transparency_thresh = -1;
    slc.fwd_model.boundary = slc.fwd_model.elems;
    show_fem(slc,[0 1 0]);
    slc = mdl_slice_mesher(img, [inf inf 2.5]);
    slc.calc_colours.transparency_thresh = -1;
    slc.fwd_model.boundary = slc.fwd_model.elems;
    show_fem(slc);
    slc = mdl_slice_mesher(img, [inf inf 1.3]);
    slc.calc_colours.transparency_thresh = -1;
    slc.fwd_model.boundary = slc.fwd_model.elems;
    show_fem(slc);
    zlim([0 3]);
    view(3)
    hold off
    subplot(133)
    hold on
    mdl_slice_mesher(img, [0 inf inf]);
    mdl_slice_mesher(img, [inf inf 2]);
    mdl_slice_mesher(img, [inf inf 2.5]);
    mdl_slice_mesher(img, [inf inf 1.3]);
    axis equal
    axis tight
    view(3)
    
    % test multi-column image
    img.elem_data(:,2) = img.elem_data;
    slc = mdl_slice_mesher(img, [0 inf inf]);
