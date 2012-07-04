function mdl = mdl_slice_mesher(fmdl,level)

if ischar(fmdl) && strcmp(fmdl,'UNIT_TEST'); do_unit_test; return, end;

mdl = fmdl;
[edges edge2elem] = get_edges(mdl);
[nodeval nodedist] = nodes_above_or_below(mdl,level);
idx = sum(nodeval(edges),2)==0; % crossed edges
crsd_edg = edges(idx,:);

dist = (nodedist(edges(idx,2)) - nodedist(edges(idx,1)));
t = -nodedist(edges(idx,1))./dist;
nodes = mdl.nodes(edges(idx,1),:) + ...
    repmat(t,1,3).*(mdl.nodes(edges(idx,2),:) - mdl.nodes(edges(idx,1),:));

% mld.nodes = level_model(fmdl, level);
show_fem(mdl);
hold on
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'bo','MarkerFaceColor','b')
hold off

function [nodeval dist] = nodes_above_or_below(mdl,level)
% nodeval = zeros(size(mdl.nodes,1),1);
% nodeval = nodeval + ((mdl.nodes(:,3) - eps) > level);
% nodeval = nodeval - ((mdl.nodes(:,3) + eps) < level);
dist = mdl.nodes(:,3) - level;
nodeval = sign(dist);

function [edges edge2elem] = get_edges(mdl)
opt.edges = true;
opt.edge2elem = true;
mdl = fix_model(mdl, opt);
edges = mdl.edges;
edge2elem = mdl.edge2elem;


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
    level = [inf 0 inf];
    mdl_slice_mesher(imdl.fwd_model, 1.3);