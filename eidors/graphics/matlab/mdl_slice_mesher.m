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
%   LEVEL  - any single-slice specification accepted by LEVEL_MODEL_SLICE
%
% Additional options can be specified as fields of MDL3D.mdl_slice_mesher
% or MDL3D.fwd_model.mdl_slice_mesher:
%       .interp_elems (defualt: true) : average element values for faces
%                                       co-planar with LEVEL
%
% To control the transparency use transparency_tresh (see CALC_COLOURS for
% details), e.g.:
%    img2d.calc_colours.transparency_thresh = -1; (no transperency)
%    calc_colours('transparency_thresh', 0.25); (some transparency)
%
% See also: LEVEL_MODEL_SLICE, SHOW_FEM, MDL_SLICE_MAPPER, SHOW_3D_SLICES, 
%           CROP_MODEL, CALC_COLOURS, PATCH

% (C) 2012-2021 Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

% TODO: 
%  1. More intuitive cut plane specification


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

fmdl.mdl_slice_mesher = parse_opt(fmdl);

opt.cache_obj = {fmdl.nodes, fmdl.elems, fmdl.mdl_slice_mesher, level};
if isfield(fmdl,'electrode');
    opt.cache_obj(end+1) = {fmdl.electrode};
end
opt.fstr      = 'mdl_slice_mesher';
opt.log_level = 4;

[nmdl, f2c, n2n, p_struct] = eidors_cache(@do_mdl_slice_mesher,{fmdl, level},opt);
nimg = build_image(nmdl, f2c, n2n, img);

switch nargout
   case 2
      out = draw_patch(p_struct, nimg.fwd_model, get_img_data(img), varargin{:});
   case 0
      out = draw_patch(p_struct, nimg.fwd_model, get_img_data(img), varargin{:});
      cmap_type = calc_colours('cmap_type');
      try 
         calc_colours('cmap_type',varargin{1}.calc_colours.cmap_type);
      end
      colormap(calc_colours('colourmap'));
      patch(out);
      calc_colours('cmap_type',cmap_type);
      clear nimg;
end


function opt = parse_opt(fmdl)
opt.interp_elems = true;

if isfield(fmdl,'mdl_slice_mesher')
    flds = fieldnames(fmdl.mdl_slice_mesher);
    for i = 1:numel(flds)
        opt.(flds{i}) = fmdl.mdl_slice_mesher.(flds{i});
    end
end


% This is a start  of a function to extrude_3d_if_reqd, so that
%   we can show 3d shapes. Currently not working
% fmdl = extrude_3d_if_reqd( fmdl );
function fmdl = extrude_3d_if_reqd( fmdl );
   if size(fmdl.nodes,2)==3; return; end
   nn = fmdl.nodes; N= size(nn,1);
   ee = fmdl.elems; E= size(ee,1);
   oN= ones(N,1);
   oE= ones(E,1) + N;
   fmdl.nodes = [nn,-10*oN; %10 is arbitrary == a guess
                 nn,+10*oN]; 
   fmdl.elems = [ee(:,[1,2,3,1]) + oE*[0,0,0,1];
                 ee(:,[1,2,3,2]) + oE*[1,0,0,1];
                 ee(:,[1,2,3,3]) + oE*[1,1,0,1]];

function [nmdl, f2c, n2n, out] = do_mdl_slice_mesher(fmdl,level)

mdl = fmdl;
opt.edge2elem = true;
opt.node2elem = true;
mdl = fix_model(mdl,opt);
edges = mdl.edges;
edge2elem = mdl.edge2elem;
tmp = mdl; 
tmp = level_model_slice( tmp, level )';
[nodeval nodedist] = nodes_above_or_below(tmp,0);
% find which edges are on electrodes
e_nodes = zeros(length(mdl.nodes),1);
try
   for i = 1:length(mdl.electrode)
      e_nodes(mdl.electrode(i).nodes) = i;
      if isfield(mdl.electrode(i),'faces')
         e_nodes(unique(mdl.electrode(i).faces)) = i;
      end
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

n2n = sparse(edges(idx,:),uint32((1:length(t))' * ones(1,2)),[1-t,t],size(mdl.nodes,1),size(nodes,1));

% nn indexes the just-created nodes, els indexes elements
if any(idx)
    [nn els] = find(edge2elem(idx,:));
else
    nn = []; els = [];
end
els_edge = els;

electrode_node = e_edges(idx);
%% crossed nodes
idx = find(nodeval == 0);
ln = length(nodes); %store the size
nodes = [nodes; mdl.nodes(idx,:)]; % add the crossed nodes to the new model
n2n = [n2n, sparse(idx,(1:length(idx)),1,size(mdl.nodes,1),length(idx))];
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

if n_tri == 0 
    error('EIDORS:NoIntersection',... 
        'No intersection found between the cut plane [%.2f %.2f %.2f] and the model.', ...
        level(1),level(2),level(3));
end
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
[nmdl.elems, idx] = sortrows(nmdl.elems);
f2c = f2c(:,idx);
[nmdl.elems, n, idx] = unique(nmdl.elems, 'rows');
if ~fmdl.mdl_slice_mesher.interp_elems
    f2c = f2c(:,n);
else
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
end

% add electrodes
try
   for i = 1:length(mdl.electrode)
      nmdl.electrode(i) = mdl.electrode(i);
      nmdl.electrode(i).nodes = find(electrode_node == i);
   end
end
out.elem_map = false(size(mdl.elems,1),1);
out.elem_map(uels) = true;
out.node_map = n2n;
out.els  = els;
out.nn   = nn;


function nimg = build_image(nmdl, f2c, n2n, img)
nimg = mk_image(nmdl,1);
if isempty(nimg.elem_data) % plane doesn't cut model
    return
end

img_data = get_img_data(img);
if size(img_data,1) == num_nodes(img.fwd_model)
    nimg.node_data = n2n' * img_data;
    nimg = rmfield(nimg, 'elem_data');
else
    nimg.elem_data = f2c' * img_data;
end
try
   nimg.calc_colours = img.calc_colours;
end


function out = draw_patch(in, nmdl, img_data, varargin)


els  = in.els;
nn   = in.nn;
nodes = nmdl.nodes;
uels = find(in.elem_map);

out.Vertices = nodes;
out.Faces    = NaN(length(uels),4);

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
    out.Faces(i,1:length(id)) = nn(id);   
end
switch length(img_data)
    case size(in.elem_map,1)
        data = img_data(uels);
        out.FaceColor = 'flat';
    case size(in.node_map,1)
        data = in.node_map' * img_data ;
        out.FaceColor = 'interp';
    otherwise
        error('wrong size of data');
end
[out.FaceVertexCData, scl_data] = calc_colours(data,varargin{:});
try
   out.FaceVertexAlphaData = double(abs(scl_data) > varargin{1}.calc_colours.transparency_thresh);
   out.FaceAlpha = 'flat';
end

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
    mdl_slice_mesher(img, struct('centre',[-.5,-.5,2],'normal_angle',[1 1 .2]))
    axis equal
    axis tight
    view(3)
    
    % test multi-column image
    img.elem_data(:,2) = img.elem_data;
    slc = mdl_slice_mesher(img, [0 inf inf]);
