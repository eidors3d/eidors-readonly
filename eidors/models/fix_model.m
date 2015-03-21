function [mdl] = fix_model(mdl,opt)
% FIX_MODEL: Add useful fields to a model
%    [mdl] = fix_model(mdl,options)
% INPUT:
%    mdl - an FEM model with at least the following fields:
%       .name
%       .nodes
%       .elem
%    options - a struct with logical values specifying which fields to
%        compute. Defaults to false for absent fields. 
%
% Run fix_model('options') to get an all-false options struct, or
% fix_model('options',true) to get an all-true one.
% mdl.faces is only replaced if necessary. mdl.boundary is never replaced.
%
% OUTPUT:
%    mdl - a copy of the input model with these additional fields:
%       .boundary
%       .boundary_face
%       .faces
%       .face2elem
%       .elem2face
%       .elem_centre
%       .face_centre
%       .face_area
%       .normals
%       .max_edge_len (per elem)
%       .edges
%       .edge_length
%       .edge2elem
%       .elem2edge
%       .face2edge
%       .node2elem
%
% For triangular meshes, edges and faces are the same.
%
% The elems will be reordered so that all faces are counter-clockwise.

% (C) 2011 Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$

if ischar(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end

if ischar(mdl) && strcmp(mdl,'options'); 
   if nargin < 2, opt = false; end
   mdl = list_options(opt); 
   return 
end
doall = false;
if nargin > 1
   opt = fix_options(mdl,opt);
else
   doall = true;
end


if doall || opt.boundary
   if ~isfield(mdl,'boundary')
      mdl.boundary = find_boundary(mdl);
   end
end
if doall || opt.faces
   if doall || opt.linear_reorder
      mdl = linear_reorder(mdl); %counter-clockwise
   end
   if elem_dim(mdl) == mdl_dim(mdl);
      [mdl.faces mdl.elem2face] = calc_faces(mdl);
   else
      % surface mesh
      mdl.faces = sort(mdl.elems,2);
      mdl.elem2face = (1:length(mdl.faces))';
   end
end
if doall || opt.face2elem
    mdl.face2elem = calc_face2elem(mdl.elem2face);
end
if doall || opt.boundary_face
   mdl.boundary_face = mdl.face2elem(:,2)==0;
end
if doall || opt.elem_centre
   mdl.elem_centre = interp_mesh(mdl, 0);
end
if doall || opt.face_centre
   tmp = mdl;
   tmp.elems = tmp.faces;
   mdl.face_centre = interp_mesh(tmp,0);
end
if doall || opt.normals
   mdl.normals = calc_normals(mdl);
end
if doall || opt.inner_normal
   if mdl_dim(mdl) == elem_dim(mdl);
      mdl.inner_normal = test_inner_normal( mdl );
   else
      eidors_msg('@@@ Inner normal test for surface meshes not implemented.',1);
   end
end
if doall || opt.max_edge_len
   mdl.max_edge_len = calc_longest_edge(mdl.elems,mdl.nodes);
end
if doall || opt.elem_volume
   mdl.elem_volume = get_elem_volume(mdl);
end
if doall || opt.face_area
   if elem_dim(mdl) == 2
      mdl.face_area = mdl.elem_volume;
   else
      mdl.face_area = calc_face_area(mdl);
   end
end
el_dim = elem_dim(mdl);
if doall || opt.edges
    if mdl_dim(mdl)==3 %that's unlikely to work for higher order elements
        [mdl.edges mdl.elem2edge] = calc_faces(mdl,2);
    else 
        mdl.edges = mdl.faces;
        mdl.elem2edge = mdl.elem2face;
    end
end
if doall || opt.node2elem
    mdl.node2elem = calc_edge2elem(mdl.elems,size(mdl.nodes,1));
end
if doall || opt.edge2elem
    mdl.edge2elem = calc_edge2elem(mdl.elem2edge, size(mdl.edges,1));
end

if doall || opt.face2edge
   if el_dim <3
      mdl.face2edge = mdl.elem2edge;
   else
      mdl.face2edge = uint32(calc_face2edge(mdl));
   end
end

if doall || opt.edge_length
   mdl.edge_length = calc_edge_length(mdl);
end
   

% decrease memory footprint
mdl.elems = uint32(mdl.elems);
if doall || opt.faces
   mdl.faces = uint32(mdl.faces);
   mdl.elem2face = uint32(mdl.elem2face);
end
if doall || opt.face2elem
   mdl.face2elem = uint32(mdl.face2elem);
end     
if doall || opt.edges
   mdl.edges = uint32(mdl.edges);
   mdl.elem2edge = uint32(mdl.elem2edge);
end
if doall || opt.edge2elem
   mdl.edge2elem = logical(mdl.edge2elem);
end     

% Test whether normal points into or outside
% mdl.inner_normal(i,j) = 1 if face j of elem i points in
function inner_normal = test_inner_normal( mdl )
   inner_normal = false(size(mdl.elem2face));
   d = elem_dim(mdl) + 1;
   for i=1:num_elems(mdl);
      el_faces = mdl.elem2face(i,:);
      el_ctr   = repmat( mdl.elem_centre(i,:), d, 1);
      vec_fa_el= el_ctr -  mdl.face_centre(el_faces,:);
      normal_i  = mdl.normals(el_faces,:);
      dot_prod = sum( normal_i.*vec_fa_el, 2 );
      inner_normal(i,:) = dot_prod' > 0;
   end

function [faces,  elem2face] = calc_faces(mdl, facedim)

e_dim = elem_dim(mdl);
if nargin == 1
    facedim = e_dim;
end

idx = nchoosek(1:e_dim+1, facedim);
elem_sorted = sort(mdl.elems,2);
[faces ib ia] = unique(reshape(elem_sorted(:,idx),[],facedim),'rows');
elem2face = reshape(ia,[],size(idx,1));

function edge2elem = calc_edge2elem(elem2edge,n_edges)

    [n_elems, el_faces] = size(elem2edge);
    elem2edgeno = (1:n_elems)'*ones(1,el_faces);
    elem2edgeno = elem2edgeno(:);
    elem2edge   = elem2edge(:);
    edge2elem = sparse(elem2edge,elem2edgeno,ones(size(elem2edgeno)),n_edges,n_elems);

function f2e = calc_face2edge(mdl)
%faces and edges are both row-wise sorted
nf = length(mdl.faces);
list(1:3:3*nf,:) = mdl.faces(:,1:2);
list(2:3:3*nf,:) = mdl.faces(:,[1 3]);
list(3:3:3*nf,:) = mdl.faces(:,2:3);
[jnk f2e] = ismember(list, mdl.edges, 'rows');
f2e = reshape(f2e,3,[])';

function face2elem = calc_face2elem(elem2face)
% This is easier to understand but very slow
%     n_face = max(elem2face(:));
%     face2elem = zeros(n_face,2);
%     for i= 1:n_face
%         [el jnk] = find(elem2face==i);
%         if numel(el)==1, el(2) = 0; end
%         face2elem(i,:) = el;
%     end
%     bck = face2elem; face2elem=[];
    [n_elems, el_faces] = size(elem2face);
    elem2faceno = (1:n_elems)'*ones(1,el_faces);
    elem2faceno = elem2faceno(:);
    elem2face   = elem2face(:);
    face2elem(elem2face,2) = elem2faceno;
    % flipping will give us the other element for shared faces
    elem2faceno = flipud(elem2faceno);
    elem2face   = flipud(elem2face);
    face2elem(elem2face,1) = elem2faceno;
    % replace with zeros repeated entries (boundary faces)
    face2elem( face2elem(:,1) == face2elem(:,2), 2) = 0;

% This function is obsolete 
function elem2face = calc_elem2face(face2elem, face_per_elem)
    n_elem = max(face2elem(:));
    elem2face = zeros(n_elem,face_per_elem);
    for i = 1:n_elem
        [f jnk] = find(face2elem==i);
        elem2face(i,:) = f;
    end
    
function normals = calc_normals(mdl)
    [n_faces face_dim] = size(mdl.faces);
    switch face_dim
        case 2
            A = mdl.nodes(mdl.faces(:,1),:);
            B = mdl.nodes(mdl.faces(:,2),:);
            normals = (B-A)*[0 1; -1 0];
        case 3
            % vectorise cross product
            x1 = mdl.nodes(mdl.faces(:,2),1) - mdl.nodes(mdl.faces(:,1),1);
            y1 = mdl.nodes(mdl.faces(:,2),2) - mdl.nodes(mdl.faces(:,1),2);
            z1 = mdl.nodes(mdl.faces(:,2),3) - mdl.nodes(mdl.faces(:,1),3);
            x2 = mdl.nodes(mdl.faces(:,3),1) - mdl.nodes(mdl.faces(:,1),1);
            y2 = mdl.nodes(mdl.faces(:,3),2) - mdl.nodes(mdl.faces(:,1),2);
            z2 = mdl.nodes(mdl.faces(:,3),3) - mdl.nodes(mdl.faces(:,1),3);
            %(a2b3 ? a3b2, a3b1 ? a1b3, a1b2 ? a2b1).
            normals = zeros(n_faces,3);
            normals(:,1) = y1.*z2 - z1.*y2;
            normals(:,2) = z1.*x2 - x1.*z2;
            normals(:,3) = x1.*y2 - y1.*x2;
        otherwise;
            error('not 2D or 3D')
    end
    normals = normals./ repmat(sqrt(sum(normals.^2,2))',face_dim,1)';
    
 function A = calc_face_area(mdl)
A = mdl.nodes(mdl.faces(:,2),:) - mdl.nodes(mdl.faces(:,1),:);
B = mdl.nodes(mdl.faces(:,3),:) - mdl.nodes(mdl.faces(:,1),:);
A = sqrt(sum(cross3(A,B).^2,2))/2;

function L = calc_edge_length(mdl)
L = sqrt(sum( (mdl.nodes(mdl.edges(:,1),:) ...
             - mdl.nodes(mdl.edges(:,2),:) ).^2 ,2 ));
    
function len = calc_longest_edge(elems,nodes)
    [E_num E_dim] = size(elems);

    pairs = nchoosek(1:E_dim,2);
    len = zeros(E_num,1);
    for i = 1:size(pairs,1)
        a = nodes(elems(:,pairs(i,1)),:);
        b = nodes(elems(:,pairs(i,2)),:);
        tmp = sqrt(sum((a-b).^2,2));
        len = max(len,tmp);  
    end
    
function out = fix_options(mdl, opt)
    out = list_options(false);
    out.linear_reorder = 1;
    flds = fieldnames(opt);
    for i = 1:length(flds)
       try
       out.(flds{i}) = opt.(flds{i});
       catch
          warning(sprintf('Option %s not recognised. Ignoring', flds{i}));
       end
    end 
    if out.boundary_face
       out.face2elem = true;
    end
    if out.inner_normal
       out.normals = true;
       out.elem_centre = true;
       out.face_centre = true;
    end
    if any([ out.boundary_face out.face_centre out.normals]) && ~isfield(mdl,'faces')
          out.faces = true;
    end
    if any([out.face2elem out.elem2face out.face_area])
       out.faces = true;
    end
    if any([out.edge2elem out.elem2edge out.edge_length out.face2edge])
        out.edges = true;
    end
    if out.edges && elem_dim(mdl) < 4
        out.faces = true;
    end

    
function out = list_options(val)
    nodes = [0 0; 0 1; 1 1; 1 0];
    elems = [1 2 3; 1 3 4];
    mdl = eidors_obj('fwd_model','square','nodes', nodes, 'elems', elems);
    out = fix_model(mdl);
    out = rmfield(out,{'elems','nodes','name','type'});
    flds = fieldnames(out);
    for i = 1:length(flds)
       out.(flds{i}) = val;
    end
    
function do_unit_test
    % square
    nodes = [0 0; 0 1; 1 1; 1 0];
    elems = [1 2 3; 1 3 4];
    mdl = eidors_obj('fwd_model','square','nodes', nodes, 'elems', elems);
    out = fix_model(mdl);
    unit_test_cmp('2d: faces'    ,out.faces    ,[1,2;1,3;1,4;2,3;3,4]);
    unit_test_cmp('2d: elem2face',out.elem2face,[1,2,4;2,3,5]);
    unit_test_cmp('2d: face2elem',out.face2elem,[1,0;2,1;2,0;1,0;2,0]);

    %cube
    nodes = [0 0 0; 0 1 0; 1 1 0; 1 0 0;...
             0 0 1; 0 1 1; 1 1 1; 1 0 1];
    elems = [1 2 3 6; 3 6 7 8; 1 5 6 8; 1 3 4 8; 1 3 6 8];
    mdl = eidors_obj('fwd_model','cube','nodes', nodes, 'elems', elems);
    out = fix_model(mdl);
    unit_test_cmp('3d: faces'    ,out.faces(1:4,:), [1,2,3;1,2,6;1,3,4;1,3,6]);
    unit_test_cmp('3d: elem2face',out.elem2face, [1,2,4,10; 12,13,14,16;
            7,8,9,15; 3,5,6,11; 4,5,9,13]);
    unit_test_cmp('3d: face2elem',out.face2elem(1:5,:), [1,0; 1,0; 4,0; 5,1; 4,5]);

    
    % test options
    opt = fix_model('options',false);
    flds = fieldnames(opt);
    % if there are no errors, option interdependence is dealt with
    % correctly
    for i = 1:length(flds)
       opt.(flds{i}) = true;
       out = fix_model(mdl, opt);
       opt.(flds{i}) = false;
    end
    
    mdl = mk_common_model('n3r2',[16,2]); mdl= mdl.fwd_model;
    out = fix_model(mdl);
    
