function [mdl] = fix_model(mdl,options)
% FIX_MODEL: Add useful fields to a model
%    [mdl] = fix_model(mdl,options)
% INPUT:
%    mdl - an FEM model with at least the following fields:
%       .name
%       .nodes
%       .elem
%    options - not coded yet
%
% OUTPUT:
%    mdl - a copy of the input model with these additional fields:
%       .boundary
%       .boundary_face
%       .faces
%       .face2elem
%       .elem2face
%       .elem_centre
%
% The elems will be reordered so that all faces are counter-clockwise.

% (C) 2011 Bartlomiej Grychtol. Licensed under GPL v2 or v3
% $Id$

if isstr(mdl) && strcmp(mdl,'UNIT_TEST'); do_unit_test; return; end


mdl = linear_reorder(mdl); %counter-clockwise
if ~isfield(mdl,'boundary')
    mdl.boundary = find_boundary(mdl);
end
[mdl.faces mdl.face2elem mdl.elem2face] = calc_faces(mdl);
mdl.boundary_face = mdl.face2elem(:,2)==0;
mdl.elem_centre = interp_mesh(mdl, 0);
tmp = mdl;
tmp.elems = tmp.faces;
%mdl.face_centre = interp_mesh(tmp,0);
mdl.normals = calc_normals(mdl);





% decrease memory footprint
mdl.elems = uint32(mdl.elems);
mdl.faces = uint32(mdl.faces);
mdl.elem2face = uint32(mdl.elem2face);
mdl.face2elem = uint32(mdl.face2elem);


function [faces face2elem elem2face] = calc_faces(mdl)

faces = [];
n_elem = num_elems(mdl);
e_dim = mdl_dim(mdl);
f_dim = e_dim-1;
n_face = e_dim+1; 

%elem2face = zeros(size(mdl.elems));
switch e_dim
    case 2
        idx = [1 2; 2 3; 1 3];
    case 3
        idx = [1 2 3; 2 3 4; 1 2 4; 1 3 4];
end
elem_sorted = sort(mdl.elems,2);
[faces ib ia] = unique(reshape(elem_sorted(:,idx),[],e_dim),'rows');
elem2face = reshape(ia,[],e_dim+1);
face2elem = calc_face2elem(elem2face);

function face2elem = calc_face2elem(elem2face)
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
    elem2faceno = flipud(elem2faceno);
    elem2face   = flipud(elem2face);
    face2elem(elem2face,1) = elem2faceno;
    face2elem( face2elem(:,1) == face2elem(:,2), 2) = 0;

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
    end
    normals = normals./ repmat(sqrt(sum(normals.^2,2))',face_dim,[])';
    
    
    
    
    
function do_unit_test
    % square
    nodes = [0 0; 0 1; 1 1; 1 0];
    elems = [1 2 3; 1 3 4];
    mdl = eidors_obj('fwd_model','square','nodes', nodes, 'elems', elems);
    out = fix_model(mdl);
    out.faces
    out.elem2face
    out.face2elem

    %cube
    nodes = [0 0 0; 0 1 0; 1 1 0; 1 0 0;...
             0 0 1; 0 1 1; 1 1 1; 1 0 1];
    elems = [1 2 3 6; 3 6 7 8; 1 5 6 8; 1 3 4 8; 1 3 6 8];
    mdl = eidors_obj('fwd_model','cube','nodes', nodes, 'elems', elems);
    out = fix_model(mdl);
    out.faces
    out.elem2face
    out.face2elem    
    
    mdl = mk_common_model('n3r2',16); mdl= mdl.fwd_model;
    out = fix_model(mdl);
    