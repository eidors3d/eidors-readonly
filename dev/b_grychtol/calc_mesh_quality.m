function Q = calc_mesh_quality(mdl)
%CALC_MESH_QUALITY Various measures of mesh quality.
%  Q = CALC_MESH_QUALITY(MDL) calculates several measures of mesh quality
%  for the EIDORS fwd_model struct MDL and returns them in a struct Q:
%     Q.NSR    - inradius to circumradius ratio (normalized shape ratio)
%                 n*r / R
%                where r = inradius; R = circumradius; n = dimensions [3,4]
%     Q.mu     - inradius to longest edge ratio (L)                     [5]
%                 r / L
%     Q.tau    - shortest to longest edge ratio                         [5]
%                 l/L
%     Q.reg    - "regularity" of tetrahedron                      [6,p.151]
%                 4r/H
%                where H is the longest altitude
%     Q.zeta   - measure relating volume and area                       [7]
%                 V^4 / ( SUM( Ai^2 ) ^3 )
%                where V = volume; Ai = area of face i
%     Q.eta    - measure realting volume and edge length                [7]
%                 V^(2/3) / SUM( li^2 )
%                where li = length of edge i
%     Q.alpha  - ratio of volume to mean edge length                    [6]
%                 V / (1/6 SUM( Ai ))
%                where Ai = area of face i
%     Q.gamma  - ratio of volume to RMS edge length                     [8]
%                 V / SQRT( 1/6 SUM( Ai^2 ) )
%
%  All the above measures are scaled to be bounded between 0 (worst) and 
%  1 (best), where a regular tetrahedron scores 1, and non-dimensional.
%
%  Measures where implemented based on review articles of V. Phathasarahy
%  (1994) [1] and D. Field (2000) [2]. The references above indicate the
%  original sources, although those were not consulted.
%
%  REFERENCES:
%  [1] Parthasarathy V (1994) "A comparison of tetrahedron quality measures" 
%      Finite Elem Anal Des 15:255–261
%  [2] Field D (2000) "Qualitative measures for initial meshes" Int J Numer
%      Meth Eng, 906:887–906
%  [3] Cavendish JC, Field DA, Frey WH (1985) "An approach to automatic
%      three-dimensional finite element mesh generation" Int J Numer Met 
%      Eng 21:329-47
%  [4] Field DA (1991) "A generic Delaunay triangulation algorithm for
%      finite elemen meshes" Adv Eng Softw 13:263-72
%  [5] Baker TJ (1989) "Element quality in tetrahedral meshes" Proc 7th Int
%      Conf on Finite Element Methods in Flow Problems, Huntsville, AL.
%      1018-24
%  [6] Dannelongue HH, Tanguy PA (1991) "Three-dimensional adaptive finite 
%      element computations and applications to non-newtonian flows" Int J 
%      Numer Meth Fl 13:145–165
%  [7] Cougny HL, Shephard MS, Georges MK (1990) "Explicit node point
%      smoothing within Octree" Report No. 10-1990, SCOREC, RPI, Troy, NY.
%  [8] Parthasarathy V, Kodiyalam (1991) "A constrained optimization
%      approach to finite element mesh smoothing" Finite Elem Anal Des
%      9:309-20
%
% See also FIX_MODEL

% TODO: Add surface mesh measures

% (C) 2013 Bartlomiej Grychtol. License: GPL v2 or v3.
% $Id$


if ischar(mdl) && strcmp(mdl, 'UNIT_TEST'), do_unit_test; return; end

opt.elem2edge = 1;
opt.face2elem = 1;
opt.boundary  = 1;
opt.normals   = 1;
opt.face_centre = 1;

mdl = fix_model(mdl, opt);
mdl.edge_length = edge_length(mdl);
mdl.max_edge_length = max(mdl.edge_length(mdl.elem2edge),[],2);
mdl.min_edge_length = min(mdl.edge_length(mdl.elem2edge),[],2);
mdl.elem_vol    = elem_volume(mdl);
mdl.face_area   = face_area(mdl);
mdl.elem_area   = elem_area(mdl);
mdl.rad_inscr   = inscribed_radius(mdl);
mdl.rad_circ    = circumscribed_radius(mdl);
mdl.tet_alt     = tet_altitudes(mdl);

% Quality measures
Q.NSR          = normalized_shape_ratio(mdl);
% this measure is not very sensitive and not bound between 0 and 1
% Q.omega        = omega(mdl);
Q.mu           = mu(mdl);
Q.tau          = tau(mdl);
Q.reg          = R_measure(mdl);
% These two measures from Field's paper are not scale invariant
% Q.R_star       = R_star(mdl);
% Q.S            = S_measure(mdl);
Q.zeta         = zeta(mdl);
Q.eta          = eta(mdl);
Q.alpha        = alpha(mdl);
Q.gamma        = gamma(mdl);


function G = gamma(mdl)
G = 6 * sqrt(2) * mdl.elem_vol ./ sqrt(mean(mdl.edge_length(mdl.elem2edge).^2,2)).^3 ;

function A = alpha(mdl)
A = 6 * sqrt(2) * mdl.elem_vol ./ mean(mdl.edge_length(mdl.elem2edge),2).^3 ;

function E = eta(mdl)
E = 12 * (3*mdl.elem_vol).^(2/3) ./ sum(mdl.edge_length(mdl.elem2edge).^2,2);

function Z = zeta(mdl)
Z = 3^7 * mdl.elem_vol.^4 ./ ...
    sum(mdl.face_area(mdl.elem2face).^2 ,2).^3;

function S = S_measure(mdl)
S = mdl.elem_vol ./ sum( mdl.edge_length(mdl.elem2edge).^2,2);

function R = R_star(mdl)
R = sqrt(2)*mdl.elem_vol ./ sum(mdl.edge_length(mdl.elem2edge),2);

function R = R_measure(mdl)
R = 4*mdl.rad_inscr ./ max(mdl.tet_alt,[],2);

function H = tet_altitudes(mdl)
elem_sorted = sort(mdl.elems,2);
v = [4 3 2 1]; % which vertex is not on the face
for i = 1:4 % loop over vertices
      Nf = mdl.normals(mdl.elem2face(:,i),:);
      Cf = mdl.face_centre(mdl.elem2face(:,i),:);
      Pe = mdl.nodes(elem_sorted(:,v(i)), :);
      H(:,i) = abs(my_dot(Nf, Pe - Cf));
end
function t = tau(mdl)
t = mdl.min_edge_length ./ mdl.max_edge_length;


function m = mu(mdl)
m = 2*sqrt(6)*mdl.rad_inscr./mdl.max_edge_length;

function O = omega(mdl)
O = sqrt(3)*mdl.max_edge_length ./ (2*sqrt(2)*mdl.rad_circ);

function NSR = normalized_shape_ratio(mdl)
NSR = 3 * (mdl.rad_inscr ./ mdl.rad_circ);


function R = circumscribed_radius(mdl)
ne = size(mdl.elems,1);
E = mdl.elems';
idx = 1:numel(E); idx(1:4:end) = [];
N = mdl.nodes(E(idx),:) - reshape(repmat(mdl.nodes(E(1:4:end),:)',3,[]),3,[])';
% prepare a matrix for solving a system
% A = sparse(length(N),length(N));
X = repmat((1:3*ne)',1,3);
T = reshape((1:3*ne)',3,[])';
Y(1:3:3*ne,1:3) = T;
Y(2:3:3*ne,1:3) = T;
Y(3:3:3*ne,1:3) = T;
A = sparse(X,Y,N);
B = sum(N.^2,2);
v = A\B;
V = reshape(v,3,[])';
R = 0.5 * sqrt(sum(V.^2,2));
C = mdl.nodes(E(1:4:end),:) + 0.5*V; % center


function R = inscribed_radius(mdl)
R = 3* mdl.elem_vol ./ mdl.elem_area;

function A = elem_area(mdl)
A = sum(mdl.face_area(mdl.elem2face),2);

function A = face_area(mdl)
F = mdl.faces';
A = mdl.nodes(mdl.faces(:,2),:) - mdl.nodes(mdl.faces(:,1),:);
B = mdl.nodes(mdl.faces(:,3),:) - mdl.nodes(mdl.faces(:,1),:);
A = sqrt(sum(my_cross(A,B).^2,2))/2;
A = A';

function V = elem_volume(mdl)
E = mdl.elems';
idx = 1:numel(E); idx(1:4:end) = [];
N = mdl.nodes(E(idx),:) - reshape(repmat(mdl.nodes(E(1:4:end),:)',3,[]),3,[])';
V = abs(my_det(N))/6;

function d = my_dot(a,b)
d = sum(a.*b,2);

function c = my_cross(a,b)
c = [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
     a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];

function D = my_det(a)
ln = size(a,1);
c1 = 1:3:ln;
c2 = 2:3:ln;
c3 = 3:3:ln;
D = a(c1,1).*a(c2,2).*a(c3,3) + ...
    a(c2,1).*a(c3,2).*a(c1,3) + ...
    a(c3,1).*a(c1,2).*a(c2,3) - ...
    a(c3,1).*a(c2,2).*a(c1,3) - ...
    a(c2,1).*a(c1,2).*a(c3,3) - ...
    a(c1,1).*a(c3,2).*a(c2,3);


function L = edge_length(mdl)
L = sqrt(sum( (mdl.nodes(mdl.edges(:,1),:) ...
             - mdl.nodes(mdl.edges(:,2),:) ).^2 ,2 ));

function do_unit_test

nodes = [0 0 0; 0 1 0; 1 1 0; 1 0 0;...
   0 0 1; 0 1 1; 1 1 1; 1 0 1];
elems = [1 2 3 6; 3 6 7 8; 1 5 6 8; 1 3 4 8; 1 3 6 8];
cube = eidors_obj('fwd_model','cube','nodes', nodes, 'elems', elems);
% note the the 5th element is a regular tet
Q = calc_mesh_quality(cube);

f = fields(Q);
for i = 1:length(f)
   disp(f{i});
   disp(Q.(f{i}));
end

real = mk_library_model('pig_23kg_16el');
Q = calc_mesh_quality(real);

f = fields(Q);
for i = 1:length(f)
   subplot(3,3,i)
   hist(Q.(f{i}),100);
   xlabel(f{i});
   xlim([0 1]);
end
subplot(3,3,9)
hist(get_elem_volume(real),100)
xlabel('volume');