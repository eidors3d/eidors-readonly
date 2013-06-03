function [Q mdl] = calc_mesh_quality(mdl, show)
%CALC_MESH_QUALITY Various measures of mesh quality.
%  [Q MDL] = CALC_MESH_QUALITY(MDL) calculates several measures of mesh 
%  quality for the EIDORS fwd_model struct MDL and returns them in Q 
%  (described below).
%  [Q MDL] = CALC_MESH_QUALITY(MDL, SHOW) if SHOW is TRUE displays two 
%  figures with histograms of each measure.
%  
%   Q.tet.
%     NSR    - inradius to circumradius ratio (normalized shape ratio)
%                 n*r / R
%                where r = inradius; R = circumradius; n = dimensions [3,4]
%     mu     - inradius to longest edge ratio (L)                     [5]
%                 r / L
%     tau    - shortest to longest edge ratio                         [5]
%                 l/L
%     reg    - "regularity" of tetrahedron                      [6,p.151]
%                 4r/H
%                where H is the longest altitude
%     zeta   - measure relating volume and area                       [7]
%                 V^4 / ( SUM( Ai^2 ) ^3 )
%                where V = volume; Ai = area of face i
%     eta    - measure realting volume and edge length                [7]
%                 V^(2/3) / SUM( li^2 )
%                where li = length of edge i
%     alpha  - ratio of volume to mean edge length                    [6]
%                 V / (1/6 SUM( Ai ))
%                where Ai = area of face i
%     gamma  - ratio of volume to RMS edge length                     [8]
%                 V / SQRT( 1/6 SUM( Ai^2 ) )
%     min_angle - minimum dihedral angle
%
%  Q.tri.
%     NSR   - inradius to circumradius ratio (normalized shape ratio)
%                 n*r / R
%                where r = inradius; R = circumradius; n = dimensions [3,4]
%     mu    - inradius to longest edge ratio (L)                        [9]
%                 r / L       = tri_mu(mdl);
%     eta   - measure relating triangle area and edge length        [10,11]
%                 A / SUM( li^2) 
%     theta - measure relating shortest altitude and edge length        [2]
%                 h / SQRT( SUM ( li^2 ) )
%     iota  - measure relating triangle area and edge length           [12]
%                 A / SUM( li ) ^2
%     kappa - ratio of shortest altitude to longest edge             [2,13]
%                 h / L
%     min_angle - minimum angle
%
%  All the above measures are non-dimensional and scaled to be bounded 
%  between 0 (worst) and 1 (best), where a regular tetrahedron/triangle
%  scores 1.
%
%  Measures where implemented based on review articles of V. Phathasarahy
%  (1994) [1] and D. Field (2000) [2]. The references above indicate the
%  original sources, although those were not consulted.
%
%  Additionally, CALC_MESH_QUALITY returns the original MDL with several
%  useful fields added.
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
%  [9] Perronnet A (1992) "Triangulation par arbre-4 de triangles equilateraux 
%      at maximisation de la qualite" Publication du Laboratoire d’Analyse
%      Numerique, Universite Pierr it Marie Curie et Centre National de la
%      Recherche Scientifique
% [10] Bhatia RP, Lawrence KL (1990) "Two-dimensional finite element mesh 
%      generation based on stripwise automatic triangulation" Comput Struct
%      36:309–319.
% [11] Bank RE, Xu J (1996) "An algorithm for coarsening unstructured 
%      meshes" Numer Math 73:1–36
% [12] Watabayshi GY, Galt JA (1986) "An optimized triangular mesh system 
%      from random points" in "Numerical Grid Generation in Computational 
%      Fluid Dynamics",Hauser J, Taylor C. (eds). Pineridge Press: Swansea,
%      1986; 437–438
% [13] Suhara J, Fukuda J (1972) "Automatic mesh generation for finite 
%      element analysis" on "Advances in Computational Methods in 
%      Structural Mechanics and Design" Oden JT, Clough RW, Yamamoto Y (eds).
%      University of Alabamal Press: Tuscaloosa, AL, 1972; 607–624
%
% See also FIX_MODEL

% (C) 2013 Bartlomiej Grychtol. License: GPL v2 or v3.
% $Id$


if ischar(mdl) && strcmp(mdl, 'UNIT_TEST'), do_unit_test; return; end

if nargin < 2
   show = 0;
end

opt.elem2edge = 1;
opt.face2elem = 1;
opt.boundary  = 1;
opt.boundary_face = 1;
opt.normals   = 1;
opt.face_centre = 1;
opt.inner_normal = 1;
mdl = fix_model(mdl, opt);

% tet properties
mdl.edge_length = edge_length(mdl);
mdl.max_edge_length = max(mdl.edge_length(mdl.elem2edge),[],2);
mdl.min_edge_length = min(mdl.edge_length(mdl.elem2edge),[],2);
mdl.elem_vol    = elem_volume(mdl);
mdl.face_area   = face_area(mdl);
mdl.elem_area   = elem_area(mdl);
mdl.rad_insphere      = insphere_radius(mdl);
mdl.rad_circumsphere  = circumsphere_radius(mdl);
mdl.tet_alt     = tet_altitudes(mdl);
% Note: The measure based on the minimum solid angle (sin(1/2 min angle)
% is bound between 0 and 1, but both 0 and 1 are degenerate, so measure is
% not calculated
mdl.solid_angle    = solid_angles(mdl);
mdl.dihedral_angle = dihedral_angles(mdl);

% tri properties
mdl.face2edge        = calc_face2edge(mdl);
mdl.rad_incircle     = incircle_radius(mdl);
mdl.rad_circumcircle = circumcircle_radius(mdl);
mdl.tri_alt          = tri_altitudes(mdl);
mdl.tri_angle        = tri_angles(mdl);

% TRI Quality measures
Q.tri.NSR          = tri_normalized_shape_ratio(mdl);
Q.tri.mu           = tri_mu(mdl);
Q.tri.eta          = tri_eta(mdl);
Q.tri.theta        = tri_theta(mdl);
Q.tri.iota         = tri_iota(mdl);
Q.tri.kappa        = tri_kappa(mdl);
Q.tri.min_angle    = tri_min_angle(mdl);

% only report measures for boundary triangles
f = fields(Q.tri);
for i = 1:length(f)
   Q.tri.(f{i}) = Q.tri.(f{i})(mdl.boundary_face);
end



% TET Quality measures
Q.tet.NSR          = tet_normalized_shape_ratio(mdl);
% this measure is not very sensitive and not bound between 0 and 1
% Q.omega        = omega(mdl);
Q.tet.mu           = mu(mdl);
Q.tet.tau          = tau(mdl);
Q.tet.reg          = R_measure(mdl);
% These two measures from Field's paper are not scale invariant
% Q.R_star       = R_star(mdl);
% Q.S            = S_measure(mdl);
Q.tet.zeta         = zeta(mdl);
Q.tet.eta          = eta(mdl);
Q.tet.alpha        = alpha(mdl);
Q.tet.gamma        = gamma(mdl);
Q.tet.min_angle    = tet_min_angle(mdl);

if show
   display_figs(Q)
end

function display_figs(Q);
f = figure;  set(f,'Name','Tetrahedron quality');
f = fields(Q.tet);
for i = 1:length(f)
   subplot(3,3,i)
   hist(Q.tet.(f{i}),100);
   xlabel(strrep(f{i},'_','\_'));
   axis tight
   xlim([0 1]);
end
% subplot(3,3,9)
% hist(180*real.dihedral_angle(:)/pi,100)
% xlabel('dihedral angle');
% xlim([0 180])

f = figure; set(f,'Name','Surface triangle quality');
f = fields(Q.tri);
for i = 1:length(f)
   subplot(3,3,i)
   hist(Q.tri.(f{i}),100);
   xlabel(strrep(f{i},'_','\_'));
   axis tight
   xlim([0 1]);
end


function a = tet_min_angle(mdl)
a = min(mdl.dihedral_angle,[],2) / acos(1/3);


function a = tri_min_angle(mdl)
a = 3*min(mdl.tri_angle,[],2)/pi;

function k = tri_kappa(mdl)
k = 2/sqrt(3) * min(mdl.tri_alt,[],2) ...
               ./ max(mdl.edge_length(mdl.face2edge),[],2);

function i = tri_iota(mdl)
i = 12*sqrt(3)*mdl.face_area ./ sum(mdl.edge_length(mdl.face2edge),2).^2;

function G = gamma(mdl)
G = 6 * sqrt(2) * mdl.elem_vol ./ sqrt(mean(mdl.edge_length(mdl.elem2edge).^2,2)).^3 ;

function A = alpha(mdl)
A = 6 * sqrt(2) * mdl.elem_vol ./ mean(mdl.edge_length(mdl.elem2edge),2).^3 ;

function E = eta(mdl)
E = 12 * (3*mdl.elem_vol).^(2/3) ./ sum(mdl.edge_length(mdl.elem2edge).^2,2);

function E = tri_eta(mdl)
E = 4 * sqrt(3) * mdl.face_area ./ sum(mdl.edge_length(mdl.face2edge).^2,2);

function Z = zeta(mdl)
Z = 3^7 * mdl.elem_vol.^4 ./ ...
    sum(mdl.face_area(mdl.elem2face).^2 ,2).^3;

function S = S_measure(mdl)
S = mdl.elem_vol ./ sum( mdl.edge_length(mdl.elem2edge).^2,2);

function R = R_star(mdl)
R = sqrt(2)*mdl.elem_vol ./ sum(mdl.edge_length(mdl.elem2edge),2);

function R = R_measure(mdl)
R = 4*mdl.rad_insphere ./ max(mdl.tet_alt,[],2);

function t = tau(mdl)
t = mdl.min_edge_length ./ mdl.max_edge_length;

function m = mu(mdl)
m = 2*sqrt(6)*mdl.rad_insphere./mdl.max_edge_length;

function m = tri_mu(mdl)
m = 2*sqrt(3)*mdl.rad_incircle./ max(mdl.edge_length(mdl.face2edge),[],2);

function t = tri_theta(mdl)
t = 2*min(mdl.tri_alt,[],2) ./ sqrt(sum(mdl.edge_length(mdl.face2edge).^2,2));

function O = omega(mdl)
O = sqrt(3)*mdl.max_edge_length ./ (2*sqrt(2)*mdl.rad_circumsphere);

function NSR = tet_normalized_shape_ratio(mdl)
NSR = 3 * (mdl.rad_insphere ./ mdl.rad_circumsphere);

function NSR = tri_normalized_shape_ratio(mdl)
NSR = 2*mdl.rad_incircle./mdl.rad_circumcircle;

function A = dihedral_angles(mdl)
v = nchoosek(1:4,2); % choose 2 faces
for i = 1:6
      N1 = mdl.normals(mdl.elem2face(:,v(i,1))         , :) ...
            .* repmat(sign(mdl.inner_normal(:,v(i,1))-0.5),1,3);
      N2 = mdl.normals(mdl.elem2face(:,v(i,2)), :) ...
            .* repmat(sign(mdl.inner_normal(:,v(i,2))-0.5),1,3);
      C  = my_cross(N1, N2);
      D  = -my_dot  (N1, N2);
      A(:,i) = atan2( sqrt(sum(C.^2,2)), D);
end


function A = solid_angles(mdl)
for i = 1:4
   E = mdl.elems';
   idx = 1:numel(E); idx(i:4:end) = [];
   N = mdl.nodes(E(idx),:) - reshape(repmat(mdl.nodes(E(i:4:end),:)',3,[]),3,[])';
   nmrtr = abs(my_det(N));
   L = sqrt(sum(N.^2,2)); % length of each vector
   dnmtr = L(1:3:end).*L(2:3:end).*L(3:3:end) ...
      + my_dot(N(1:3:end,:), N(2:3:end,:)) .* L(3:3:end) ...
      + my_dot(N(1:3:end,:), N(3:3:end,:)) .* L(2:3:end) ...
      + my_dot(N(2:3:end,:), N(3:3:end,:)) .* L(1:3:end);
   
   A(:,i) = atan2( nmrtr, dnmtr );
end
idx = A<0;
A(idx) = A(idx) + pi;
A = 2*A;

function t = tri_angles(mdl)
v = 1:3;
L = mdl.edge_length(mdl.face2edge);
L2 = L.^2;
for i = 1:3
   t(:,i) = acos( (L2(:,v(1)) + L2(:,v(2)) - L2(:,v(3)) ) ...
                    ./(2 .* L(:,v(1)) .* L(:,v(2)) ) );
   v = circshift(v,[1 2]);
end


function R = circumsphere_radius(mdl)
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


function H = tet_altitudes(mdl)
elem_sorted = sort(mdl.elems,2);
v = [4 3 2 1]; % which vertex is not on the face
for i = 1:4 % loop over vertices
      Nf = mdl.normals(mdl.elem2face(:,i),:);
      Cf = mdl.face_centre(mdl.elem2face(:,i),:);
      Pe = mdl.nodes(elem_sorted(:,v(i)), :);
      H(:,i) = abs(my_dot(Nf, Pe - Cf));
end

function H = tri_altitudes(mdl)
s = sum(mdl.edge_length(mdl.face2edge),2) / 2;
S = repmat(s,1,3);
nmrtr = 2*sqrt(s .* prod(S - mdl.edge_length(mdl.face2edge),2));
H = repmat(nmrtr,1,3) ./ mdl.edge_length(mdl.face2edge) ;



function R = circumcircle_radius(mdl)
R = 0.25*prod(mdl.edge_length(mdl.face2edge),2) ./ mdl.face_area; 

function R = incircle_radius(mdl)
R = 2*mdl.face_area ./ sum(mdl.edge_length(mdl.face2edge),2);

function R = insphere_radius(mdl)
R = 3* mdl.elem_vol ./ mdl.elem_area;

function f2e = calc_face2edge(mdl)
%faces and edges are both column-wise sorted
nf = length(mdl.faces);
list(1:3:3*nf,:) = mdl.faces(:,1:2);
list(2:3:3*nf,:) = mdl.faces(:,[1 3]);
list(3:3:3*nf,:) = mdl.faces(:,2:3);
[jnk f2e] = ismember(list, mdl.edges, 'rows');
f2e = reshape(f2e,3,[])';

function A = elem_area(mdl)
A = sum(mdl.face_area(mdl.elem2face),2);

function A = face_area(mdl)
F = mdl.faces';
A = mdl.nodes(mdl.faces(:,2),:) - mdl.nodes(mdl.faces(:,1),:);
B = mdl.nodes(mdl.faces(:,3),:) - mdl.nodes(mdl.faces(:,1),:);
A = sqrt(sum(my_cross(A,B).^2,2))/2;
A = A;

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
[Q real] = calc_mesh_quality(real, 1);



