function nimg = mdl_slice_mesher(fmdl,level)

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return, end;
switch fmdl.type
    case 'image'
        img  = fmdl;
        fmdl = fmdl.fwd_model;
    case 'fwd_model'
        img  = mk_image(fmdl,1);
    otherwise
        error('Unknown object type');
end


mdl = fmdl;
opt.edge2elem = true;
opt.node2elem = true;
mdl = fix_model(mdl,opt);
edges = mdl.edges;
edge2elem = mdl.edge2elem;
tmp = mdl;
tmp.nodes = level_model( tmp, level )';
[nodeval nodedist] = nodes_above_or_below(tmp,0);
%% crossed edges
idx = sum(nodeval(edges),2) == 0 ; 
dist = (nodedist(edges(idx,2)) - nodedist(edges(idx,1)));
t = -nodedist(edges(idx,1))./dist;
nodes = mdl.nodes(edges(idx,1),:) + ...
    repmat(t,1,3).*(mdl.nodes(edges(idx,2),:) - mdl.nodes(edges(idx,1),:));
[nn els] = find(edge2elem(idx,:));


%% crossed nodes
% idx = nodeval == 0;
% [nnn eee] = find(mdl.node2elem(idx,:));
% [ueee jnk n1] = unique(eee);
% nodes_per_elem = jnk;
% nodes_per_elem(2:end) = diff(jnk);
% % if an elem has more than 2 crossed nodes, add it
% add = find(nodes_per_elem > 2);
% els = [els; ueee(add)];
% for i = 1:length(add)
%     nn = [nn; nnn(n1 == add(i))];
% end
% [els idx] = sort(els);
% nn = nn(idx);
% % for elems with less than 4 crossed edges -> add crossed nodes if needed
% [uels jnk n] = unique(els);
% nodes_per_elem = jnk;
% nodes_per_elem(2:end) = diff(jnk);
% 
% [idx ia ib] = intersect(ueee, uels);
% for i = 1:length(ia)
%     newnodes = nnn(n1==ia(i));
%     nn = [nn; newnodes];
%     els = [els; repmat(ib(i),length(newnodes),1)];
% end
% [els idx] = sort(els);
% nn = nn(idx);
[uels jnk n] = unique(els);
nodes_per_elem = jnk;
nodes_per_elem(2:end) = diff(jnk);

n_tri = length(uels) + sum(nodes_per_elem==4);

nmdl.type = 'fwd_model';
nmdl.nodes = nodes;
nmdl.elems = zeros(n_tri,3);
nimg = mk_image(nmdl,1);
c = 1;
% TODO: Speed this up
for i = 1:length(uels)
    switch nodes_per_elem(i)
        case 3
            nmdl.elems(c,:) = nn(n==i);
            nimg.elem_data(c) = img.elem_data(uels(i));
            c = c + 1;
        case 5
            nds = nn(n==i);
            nmdl.elems(c,:) = nds(1:3);
            nimg.elem_data(c) = img.elem_data(uels(i));
            nmdl.elems(c+1,:) = nds(2:4);
            nimg.elem_data(c+1) = img.elem_data(uels(i));
            c = c + 2;
    end
end
nimg.fwd_model = nmdl;

% This bit could be useful to look at the shape of the actual elements
% But the quad patches would have to be re-ordered
% show_fem(fmdl);
% for i = 1:length(uels)
%     idx = els == uels(i);
%     patch(nodes(nn(idx),1),nodes(nn(idx),2),nodes(nn(idx),3),1)
% end


function [nodeval dist] = nodes_above_or_below(mdl,level)
% nodeval = zeros(size(mdl.nodes,1),1);
% nodeval = nodeval + ((mdl.nodes(:,3) - eps) > level);
% nodeval = nodeval - ((mdl.nodes(:,3) + eps) < level);
dist = mdl.nodes(:,3) - level;
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
    slc = mdl_slice_mesher(img, [3 3 2]);
    subplot(121)
    show_fem(img);
    subplot(122)
    show_fem(slc);
    zlim([0 3]);