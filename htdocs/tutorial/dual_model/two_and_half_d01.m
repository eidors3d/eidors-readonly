% Build 2D and 3D model $Id$

if ~exist('demo_img');
   [inhomg_img, demo_img] = demo_real;
   close all;
end

% Create 2D FEM of all NODES with z=0
f_mdl = demo_img.fwd_model;
n2d = f_mdl.nodes( ...
           (f_mdl.nodes(:,3) == 0), 1:2);
e2d = delaunayn(n2d);
c_mdl = eidors_obj('fwd_model','2d','elems',e2d,'nodes',n2d);

subplot(121);
show_fem(f_mdl); title('fine (3d) model');

subplot(122);
show_fem(c_mdl); title('coarse (2d) model');
axis square

print_convert two_and_half_d01a.png '-density 75'

% Simulate data - inhomogeneous
vi= fwd_solve(inhomg_img);

% Simulate data - homogeneous
homg_img= inhomg_img; homg_img.elem_data(:) = 1;
vh= fwd_solve(homg_img);
