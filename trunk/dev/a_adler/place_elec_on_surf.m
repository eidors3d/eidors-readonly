function mdl2 = place_elec_on_surf(mdl,elec_pos, elec_spec)
% PLACE_ELEC_ON_SURF Place electrodes on the surface of a model
%  mdl = place_elec_on_surf(mdl,elec_pos, elec_spec)

% (C) Bartlomiej Grychtol and Andy Adler, 2012. Licenced under GPL v2 or v3

% This is work in progress code. Careful!

if isstr(mdl) && strcmp(mdl, 'UNIT_TEST') do_unit_test; return; end;
oldbound = [];
try 
   oldbound = mdl.boundary;
   mdl = rmfield(mdl,'boundary');
end
mdl = fix_model(mdl);
mdl.boundary = mdl.faces(mdl.boundary_face,:);
mdl.elems = mdl.boundary;
mdl = linear_reorder(mdl,1);
mdl.boundary = mdl.elems;
mdl.faces = mdl.boundary;

mdl  = linear_reorder(mdl,1);
mdl.face_centre = mdl.face_centre(mdl.boundary_face,:);
% This is a bit convoluted
to_flip = unique(nonzeros(uint32(mdl.boundary_face(mdl.elem2face) .* mdl.inner_normal) .* mdl.elem2face));
mdl.normals(to_flip,:) = -mdl.normals(to_flip,:);
mdl.normals = mdl.normals(mdl.boundary_face,:);

mdl = rmfield(mdl, 'inner_normal');
mdl = rmfield(mdl, 'boundary_face');

elecs = parse_elecs(mdl,elec_pos,elec_spec);
for i = 1:length(elecs)
   N = grow_neighbourhood(mdl,elecs(i));
   [mdl E1{i} E2{i} V{i}] = add_electrodes(mdl,N,elecs(i));
end


STL.vertices = mdl.nodes;
STL.faces    = mdl.elems;
stl_write(STL,'tmp.stl');

call_netgen('tmp.stl','tmp.vol');
fmdl=ng_mk_fwd_model('tmp.vol',[],[],[],[]);
fmdl = fix_model(fmdl);
fmdl.boundary = fmdl.faces(fmdl.boundary_face,:);

mdl = fmdl;
mdl.elems = fmdl.boundary;

for i = 1:length(elecs)
   mdl = flatten_electrode(mdl,E1{i},E2{i}, V{i});
end
mdl = linear_reorder(mdl,1);
STL.vertices = mdl.nodes;
STL.faces    = mdl.elems;
stl_write(STL,'tmp.stl');
% TODO:
% re-mesh the inside???

mdl2 = gmsh_stl2tet('tmp.stl');
mdl2.electrode = mdl.electrode;
for i = 1:length(elecs)
   enodes = mdl.nodes(mdl.electrode(i).nodes,:);
   mdl2.electrode(i).nodes = find_matching_nodes(mdl2,enodes,1e-5);
end

function mdl = flatten_electrode(mdl,inner,outer, V)
n1 = find_matching_nodes(mdl,inner, 1e-5);
n2 = find_matching_nodes(mdl,outer, 1e-5);
% remove the side nodes of the electrode
N1 = false(length(mdl.nodes),1);
N1(n1) = true;
N2 = false(length(mdl.nodes),1);
N2(n2) = true;
rm = sum(N1(mdl.elems),2)>0 & sum(N2(mdl.elems),2)>0;



f = find(sum(N2(mdl.elems),2)>1 & ~rm,1,'first');
B = find(mdl.boundary_face);
p = mdl.face_centre(B(f),:);
r = Inf;
mdl.elems(rm,:) = [];
mdl.boundary = mdl.elems;
f = f - nnz(rm(1:f));
N = grow_neighbourhood(mdl,f,p,r);

% WARNING: Here we assume the sides of the electrode are one element high!

%nodes to move
ntm = unique(mdl.elems(N,:));
mdl.nodes(ntm,:) = mdl.nodes(ntm,:) - repmat(V,length(ntm),1);
e_nodes = ntm;

%remap outer nodes to inner ones
map = 1:length(mdl.nodes);
map(n2) = n1;
mdl.elems = map(mdl.elems);
e_nodes = map(ntm);

% remove the outer nodes
m = true(length(mdl.nodes),1);
m(n2) = false;
map = zeros(size(m));
map(m) = 1:nnz(m);

mdl.nodes(n2,:) = [];
mdl.elems = map(mdl.elems);
e_nodes = map(e_nodes);

mdl.boundary = mdl.elems;
if ~isfield(mdl,'electrode')
   mdl.electrode = struct();
   l = 1;
else
   l = length(mdl.electrode);
   % because we are changing the number of nodes, we need to correct the
   % electrodes that are there alread
   for i = 1:l
      mdl.electrode(i).nodes = map(mdl.electrode(i).nodes);
   end
   l = l + 1;
end
mdl.electrode(l).nodes = double(e_nodes);
mdl.electrode(l).z_contact = 0.01;

% img = mk_image(mdl,1);
% img.elem_data(N) = 2;
show_fem(mdl);


function match = find_matching_nodes(mdl, nodes,th)
l0 = length(mdl.nodes);
match = 0 * (1:length(nodes));
for n = 1:length(nodes)
   D = mdl.nodes - repmat(nodes(n,:),l0,1);
   D = sqrt(sum(D.^2,2));
   [val p] = min(D);
   if val < th
      match(n) = p;
   end
end

% Returns a joint surface mesh and the list of nodes on the side of the
% electrode
function [joint EL1 EL2 V] = add_electrodes(mdl,N,elecs)
fc = find_face_under_elec(mdl,elecs.pos);
% N indexes the boundary, need index into faces
% fcs = find(mdl.boundary_face);
% fcs = fcs(N);
fcs = N;

jnk.type = 'fwd_model';
jnk.elems = mdl.boundary(N,:);
jnk.nodes = mdl.nodes;
jnk.boundary = jnk.elems;
img = mk_image(jnk,1);

show_fem(jnk);
hold on
plot3(elecs.points(:,1),elecs.points(:,2),elecs.points(:,3),'ro');
% plot3(node(:,1),node(:,2),node(:,3),'bs');
% plot3(mdl.nodes(nn(outer),1), mdl.nodes(nn(outer),2), mdl.nodes(nn(outer),3),'bs') 
hold off


%nodes used
[nn,I, J] = unique(mdl.faces(fcs,:));
outer = true(size(nn));
for i = 1:length(nn)
   if sum(J==i) == sum(mdl.boundary(:) == nn(i))
      outer(i) = false;
   end
end
% we want to keep the ones on the outside
keep = false(size(mdl.nodes,1),1);
keep(nn) = outer;
keep = keep(jnk.elems);

% this will not catch the situation where the element reaches from boundary
% to boundary and the electrode is in the middle (small electrode, big
% element). Fortunately, in these cases the edge will be there twice
edges = reshape(jnk.elems(:,[1 2 2 3 3 1])',2,[])';
if size(keep,2) == 1; keep = shiftdim(keep,1); end
keep = reshape(keep(:,[1 2 2 3 3 1])',2,[])';
rm = sum(keep,2)<2;
edges(rm,:) = [];

% tmp = jnk.elems'; tmp(:,rm) = [];
% edges = reshape(tmp(keep'),2,[])'; 
% % rotate those that were 101
% tmp=edges(:,1);
% flip = keep(:,2) == 0;
% edges(flip,1) = edges(flip,2);
% edges(flip,2) = tmp(flip,1);

% detect and remove double entries
rm = ismember(edges,edges(:,[2 1]),'rows');
edges(rm,:) = [];


% project all nodes of the faces in N onto the plane of the electrode
nodes = unique(mdl.faces(fcs,:));
PN = project_nodes_on_elec(mdl,elecs,nodes);
% Ne = mdl.normals(fc,:);
% Pe = elecs.pos;
% for i = 1:length(nodes)
%    P = mdl.nodes(nodes(i),:);
%    PN(i,:) = P + dot(Pe - P, Ne) * Ne;
% end

% % electrode coordinate system
[u v s] = get_face_basis(mdl,fc);
% u = mdl.normals(fc,:); % unit normal
% % vertical vector on the plane of that surface triangle
% v = [0 0 1] - dot([0 0 1],u) *u; v = v/norm(v);
% s = cross(u,v); s= s/norm(s);

% mark nodes that are very close to elecs.points for removal
rm = false(length(PN),1);
for i = 1:length(PN)
   D = repmat(PN(i,:),length(elecs.points),1) - elecs.points;
   D = sqrt(sum(D.^2,2));
   if any(D < elecs.maxh/2)
      rm(i) = true;
   end
end
% we can only delete if it's not part of the boundary
b = unique(edges(:));
rm = find(rm);
rm(ismember(nodes(rm),b)) = [];

PN(rm,:) = [];
nodes(rm) = [];

points = [PN; elecs.points];
np = size(points,1);
x = dot(points,repmat(v,np,1),2);
y = dot(points,repmat(s,np,1),2);

map(nodes) = 1:length(nodes);
edges = map(edges); %

f = length(PN) +(1:2);
C = [];
for i= 0:length(elecs.points)-2
   C = [C; i+f];
end

D = DelaunayTri([x y],[edges; C]);
% D = DelaunayTri([x y],[edges]);
els = D.Triangulation(D.inOutStatus,:);

% project all electrode points on all triangles, using the normal of the central elem
Ne = mdl.normals(fc,:);
% % make sure the normal points outside
% if mdl.inner_normal(fc)
%    Ne = -1 * Ne;
% end
for j = 1:length(elecs.points)
   Pe = elecs.points(j,:);
   for i = 1:length(fcs)
      Nf = mdl.normals(fcs(i),:);
      Cf = mdl.face_centre(fcs(i),:);
      % the plane is (X - Cf).Nf = 0
      % the line is X = Pe + tNe (through Pe perpendicular to the main elec
      % face
      % We want X that satisfies both.
      % (Pe +tNe -  Cf).Nf = 0
      % (Pe - Cf).Nf + tNe.Nf = 0
      % t = (Cf-Pe).Nf / (Ne.Nf) 
      % X = Pe + Ne * (Cf-Pe).Nf / (Ne.Nf) 
      X = Pe + Ne * dot(Cf-Pe,Nf) / dot(Ne,Nf) ;
      if point_in_triangle(X, mdl.faces(fcs(i),:), mdl.nodes)
         Proj(j,:) = X;
         FC(j) = fcs(i);
         break;
      end
   end
end

% this is just output
EL1 = Proj;

ln = length(nodes);
IN = inpolygon(x(1:ln),y(1:ln),x(ln+1:end),y(ln+1:end));

add = elecs.maxh;

nn = mdl.nodes(nodes,:) + add * repmat(IN,1,3) .* repmat(Ne,ln,1);
le = length(elecs.points);
ne = Proj + add * repmat(Ne,le,1);

%this is just output
EL2 = ne;
V = add*Ne;

% the nodes of the electrode
IN = [IN; ones(le,1)];
el_c = D.incenters;
el_c(~D.inOutStatus,:) = [];
e_el = inpolygon(el_c(:,1),el_c(:,2),x(ln+1:end),y(ln+1:end));
% e_el =  sum(IN(els),2) == 3;
els(e_el,:) = els(e_el,:) + (els(e_el,:)>ln ) .* le;


% add connecting elements
E = [];
f = ln + [ 1 le+2 le+1; le+2 1 2];
for j = 0:(le-2)
   E = [E; j+f];
end
M = ln + [le+1 2*le le; le 1 le+1];
E = [E; M];

jnk.nodes = [nn ; Proj;  ne];
jnk.elems = [ els; E];
jnk.boundary = jnk.elems;
show_fem(jnk)

% remove the patch we're replacing
big = mdl;
big.boundary(N,:) = [];
big.faces(N,:) = [];
big.normals(N,:) = [];
big.face_centre(N,:) = [];


big.elems = big.boundary;

% big = linear_reorder(big,1);

joint = merge_meshes(big,jnk,0.001);
joint.boundary = joint.elems;
joint.faces = joint.boundary;
opt.normals = true;
opt.face_centre = true;
joint = fix_model(joint,opt);

function PN = project_nodes_on_elec(mdl,elecs,nodes)
fc = find_face_under_elec(mdl,elecs.pos);
Ne = mdl.normals(fc,:);
Pe = elecs.pos;
for i = 1:length(nodes)
   P = mdl.nodes(nodes(i),:);
   PN(i,:) = P + dot(Pe - P, Ne) * Ne;
end

% OUTPUT:
%  elecs(i).pos   = [x,y,z]
%  elecs(i).shape = 'C' or 'R'
%  elecs(i).dims  = [radius] or [width,height]
%  elecs(i).maxh  = '-maxh=#' or '';
%  elecs(i).points= list of points around the perimeter
% Angles (th) are interpreted with the mean of boundary nodes as origin
function [elecs] = parse_elecs(mdl, elec_pos, elec_shape )

if size(elec_shape,2) < 3
   elec_shape(:,3) = elec_shape(:,1)/10;
end

if size(elec_pos,1) == 1
   % Parse elec_pos = [n_elecs_per_plane,(0=equal angles,1=equal dist),z_planes]
   n_elecs= elec_pos(1); % per plane
   offset = elec_pos(2) - floor(elec_pos(2));
   switch floor(elec_pos(2))
      case 0
         th = linspace(0,2*pi, n_elecs+1)'; th(end)=[];
         th = th + offset*2*pi;
         ind = th >= 2*pi;
         th(ind) = th(ind) - 2*pi;
      case 1
         error('not implemented yet');
   end
   on_elecs = ones(n_elecs, 1);
   el_th = [];
   el_z  = [];
   for i=3:length(elec_pos)
      el_th = [el_th; th];
      el_z  = [el_z ; on_elecs*elec_pos(i)];
   end
else
   el_th = elec_pos(:,1)*2*pi/360;
   el_z  = elec_pos(:,2);
end

el_th(el_th>pi) =  el_th(el_th>pi) - 2*pi;
el_th(el_th<-pi) = el_th(el_th<-pi) + 2*pi;

n_elecs= size(el_z,1);

if size(elec_shape,1) == 1
   elec_shape = ones(n_elecs,1) * elec_shape;
end

for i = 1:n_elecs
   if 1
      [fc elecs(i).pos] = find_elec_centre(mdl,el_th(i),el_z(i));
   else
      [bfc elecs(i).pos] = find_face_under_elec(mdl,elecs(i).pos);
      idx = find(mdl.boundary_face);
      fc = idx(bfc);
   end
%    elecs(i).face = fc; % this changes too often to store!
   elecs(i).dims = elec_shape(i,1:2);
   elecs(i).dims(elecs(i).dims==0) = [];
   elecs(i).maxh = elec_shape(i,3);
   
   if elec_shape(i,2) == 0
      elecs(i).shape = 'C';
      r = elec_shape(i,1);
      n = ceil(2*pi*elec_shape(i,1) / elec_shape(i,3));
      t = linspace(0,2*pi,n); t(end) = [];
      x = r*sin(t); y = r*cos(t);
   else
      elecs(i).shape = 'R';
      height = elec_shape(i,1); width = elec_shape(i,2); d = elec_shape(i,3);
      nh = ceil(height/d); nw = ceil(height/d);
      ph = linspace(-height/2,height/2,nh);
      pw = linspace(-width/2,width/2,nw);
      y = [ph, ph(end)*ones(1,nw-2), fliplr(ph), ph(1)*ones(1,nw-2)];
      x = [pw(1)*ones(1,nh-1), pw, pw(end)*ones(1,nh-2), fliplr(pw(2:end))];
      %    % we don't want real rectangles, because Netgen will merge coplanar
      %    % faces, so we create a nice superellipse instead
      n = ceil(2*(height+width)/d);
      %    t = linspace(2*pi,0,n); t(end) = [];
      %    N = 8;
      %    x = abs(cos(t)).^(2/N) * width/2  .* sign(cos(t));
      %    y = abs(sin(t)).^(2/N) * height/2 .* sign(sin(t));
      % superellipses are also bad, what about a wavy rectange?
      [pp] = fourier_fit([x; y]',20);
      t = linspace(0,1,n+1); t(end) = [];
      xy = fourier_fit(pp,t);
      x = xy(:,1)'; y = xy(:,2)';
   end
   fc = find_face_under_elec(mdl,elecs(i).pos);
   [u v s] = get_face_basis(mdl, fc);

   
   elecs(i).points = ones(size(x))' * elecs(i).pos + x'*s + y'*v;
end

function [u v s] = get_face_basis(mdl, fc)
   u = mdl.normals(fc,:); % unit normal
   % vertical vector on the plane of that surface triangle
   v = [0 0 1] - dot([0 0 1],u) *u; v = v/norm(v);
   s = cross(u,v); s= s/norm(s);

function [fc pos] = find_elec_centre(mdl, el_th,el_z)
Ctr = mean(mdl.nodes(mdl.boundary,:));
Ctr(3) = el_z;

%1. Find elements that cross the z plane
n_above = mdl.nodes(:,3) >= el_z;
sum_above = sum(n_above(mdl.boundary),2) ;
els = sum_above > 0 & sum_above < 3;

%2. Find an edge that crosses el_th
n = unique(mdl.boundary(els,:));
nn = mdl.nodes(n,1:2);
nn = nn - repmat(Ctr(:,1:2),length(nn),1);
th = cart2pol(nn(:,1),nn(:,2));
th(:,2) = 1:length(th);
th = sortrows(th);
idx = find(th(:,1) > el_th,1,'first');
if isempty(idx) || idx == 1
   n1 = n(th(1,2));
   n2 = n(th(end,2));
else
   n1 = n(th(idx-1,2));
   n2 = n(th(idx,  2));
end
% elements in els that contain these nodes (they don't need to be on the
% same element)
el = els & sum( (mdl.boundary == n1) + (mdl.boundary == n2) ,2) > 0;
el = find(el);
% fcs = find(mdl.boundary_face);
% fcs = fcs(el);

[De(1) De(2) De(3)]  = pol2cart(el_th,1, 0); 
for i = 1:length(el)
   Nf = mdl.normals(el(i),:);
   Cf = mdl.face_centre(el(i),:);
   % the plane is (X - Cf).Nf = 0
   % the line is X = Ctr + tDe (through Ctr along De
   % We want X that satisfies both.
   % (Ctr +tDe -  Cf).Nf = 0
   % (Ctr - Cf).Nf + tDe.Nf = 0
   % t =
   % X = Ctr + De * (Cf-Ctr).Nf / (De.Nf)
   X = Ctr + De * dot(Cf-Ctr,Nf) / dot(De,Nf) ;
   if point_in_triangle(X, mdl.faces(el(i),:), mdl.nodes)
      pos = X;
      fc = el(i);
      break;
   end
   
   % project the line on this element
   % check if it falls inside
end

function write_polyhedron(fid,fmdl, str)
   fprintf(fid, 'solid %s = polyhedron(',str);
   NN = fmdl.nodes;
   BB = fmdl.boundary;
   for i = 1:size(fmdl.nodes,1)
   fprintf(fid, '\n                      %4.16f, %4.16f, %4.16f;',NN(i,:));
   end
   for i = 1:size(fmdl.boundary,1)
   fprintf(fid, ';\n                      %d, %d, %d',BB(i,:));
   end
   fprintf(fid, ');\n');

function out = grow_neighbourhood(mdl, varargin)
use_elec = false;
if length(varargin) == 1
   use_elec = true;
   elecs = varargin{1};
   fc = find_face_under_elec(mdl,elecs.pos);
   p = elecs.pos;
   switch elecs.shape
      case 'R'
         r = sqrt(sum(elecs.dims.^2,2));
      case 'C'
         r = 2 * elecs.dims(1);
   end
   % here we want fc to be an index into the boundary field
%    fc = sum(mdl.boundary_face(1:fc));
else
   fc = varargin{1};
   p = varargin{2};
   r = varargin{3};
end

done = false(length(mdl.boundary),1);
todo = false(length(mdl.boundary),1);
todo(fc) = true;
bb = mdl.boundary;
vv = mdl.nodes;
% distance of each vertex to a the line perpendicular to face fc passing
% through p
dv = vv - repmat(p,length(vv),1);
nl = mdl.normals;
nl = repmat(nl(fc,:),length(vv),1);
dd = sqrt(sum( (dv - repmat(dot(dv,nl,2),1,3) .* nl).^2,2));
dim = size(bb,2);
first = true; % at first iteration, add all neighbours
if use_elec
   PN = project_nodes_on_elec(mdl,elecs,1:length(mdl.nodes));
   emin = min(elecs.points);
   emax = max(elecs.points);
   toofar = false(size(mdl.boundary,1),1);
   
   for i = 1:3
      nodes = reshape(PN(mdl.boundary,i),[],3);
      toofar =  toofar |  sum(nodes > emax(i),2) == 3 | sum(nodes < emin(i),2) == 3;
   end
end
while any(todo)
   id = find(todo,1,'first');
   done(id) = 1;
   nn = find_neighbours(id,bb);
   if use_elec
      nn = nn & ~toofar;
   elseif first
      % include all neighbours
      first = false;
   else
      % at least one node must be close enough
      nn = nn & sum(dd(bb) <= r,2) > 0;
   end
   todo = todo | nn;
   todo(done) = 0;
%    disp(sprintf('id: %d done: %d todo: %d',id, nnz(done),nnz(todo)));
%    disp(find(todo)');
%    disp(find(done)');
end
out = find(done);


function nn =  find_neighbours(fc, bb);
dim = size(bb,2);
nn = false(length(bb),1);
for i = 1:dim
   node = bb(fc,i);
   nn = nn | sum(bb == node,2) > 0;
end
nn(fc) = 0;

function [e p] = find_face_under_elec(mdl, elec_pos)
for i = 1:size(elec_pos,1)
   ee = repmat(elec_pos(i,:),length(mdl.faces),1);
   fc = mdl.face_centre;
   n  = mdl.normals;
   proj = ee - repmat(dot(ee-fc, n,2),1,3) .* n;
   in = point_in_triangle(proj,mdl.faces,mdl.nodes);
   if nnz(in) == 1
      e(i) = find(in);  % this should be an index into mdl.boundary
      p(i,:) = proj(in,:);
   else
      % take the element that is closest to ee
      cand = find(in);
      dd = sqrt(sum(  ...
         (fc(cand,:) - repmat(elec_pos(i,:),length(cand),1)).^2,2));
      [jnk pos] = min(dd);
      e(i) = cand(pos);
      p(i,:) = proj(e(i),:);
   end

end

% check if point p is in triangle E defined by indices into vertices V
function out = point_in_triangle(p, E, V)
%http://www.blackpawn.com/texts/pointinpoly/default.html
% vectors
v0 = V(E(:,3),:) - V(E(:,1),:);
v1 = V(E(:,2),:) - V(E(:,1),:);
v2 = p  - V(E(:,1),:);

% dot products
dot00 = dot(v0, v0, 2);
dot01 = dot(v0, v1, 2);
dot02 = dot(v0, v2, 2);
dot11 = dot(v1, v1, 2);
dot12 = dot(v1, v2, 2);

% barycentric coordinates
invDenom = 1 ./ (dot00 .* dot11 - dot01 .* dot01);
u = (dot11 .* dot02 - dot01 .* dot12) .* invDenom;
v = (dot00 .* dot12 - dot01 .* dot02) .* invDenom;

out = u >= 0 & v >= 0 & (u + v < 1);


function do_unit_test
xy= [ -0.89 -0.74 -0.21  0.31  0.79  0.96  0.67  0.05 -0.36 -0.97;
       0.14  0.51  0.35  0.50  0.27 -0.23 -0.86 -0.69 -0.85 -0.46]';
[fmdl] = ng_mk_extruded_model({2,xy,[4,80],},[],[]);
elec_pos = [-0.5, -0.8, 1];
% place_elec_on_surf(fmdl, elec_pos, [0.1 0 0.01]);
% place_elec_on_surf(fmdl, elec_pos, [0.15 0.1 0.01]);
place_elec_on_surf(fmdl, [16 0 1], [0.15 0.1 0.01]);
% place_elec_on_surf(fmdl, [16 0 1], [0.1 0 0.01]);