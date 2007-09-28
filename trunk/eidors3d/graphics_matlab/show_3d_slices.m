function show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% Show a 3d view of an object with many slices through it

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: show_3d_slices.m,v 1.1 2007-09-28 00:26:36 aadler Exp $

clim = [];
ref_lev= 'use_global';
rimg= calc_slices( img, z_cuts(:)*[inf,inf,1]);
c_img = calc_colours( rimg, clim, 0, ref_lev );
np= calc_colours('npoints');

xyz_max= max(img.fwd_model.nodes);
xyz_min= min(img.fwd_model.nodes);
[x,y]=meshgrid( linspace(xyz_min(1),xyz_max(1), np), ...
                linspace(xyz_min(2),xyz_max(2), np) );

ff=isnan(rimg);
c_img(isnan(rimg))= NaN;


hh=surf(x,y,c_img*0,c_img);
set(hh,'CDataMapping','direct', ... % WHY WOULD IT BE ANYTHING ELSE !!!
       'EdgeAlpha',0);

hold on
c_img(isnan(rimg))= -1e50;
[jnk,hh]=contour(x,y,c_img, [-1e49,-1e49]);
hold off

