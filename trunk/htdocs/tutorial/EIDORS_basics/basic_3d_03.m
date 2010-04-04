% Basic 3d model $Id$

clf
% Show 3D object as slices
img2.calc_colours.greylev = -0.05;
img2.calc_colours.npoints = 128;
subplot(221);
show_3d_slices(img2, [0.5,1.5,1.8,2.1]);
view(-14,13); axis tight; axis equal; zlim([0,3]);

% Show 3D object as slices
print_convert basic_3d_03a.png

subplot(221);
show_3d_slices(img2, [1,1.9], [0.5],[0.5]);
view(-14,13); axis tight; axis equal;

print_convert basic_3d_03b.png
