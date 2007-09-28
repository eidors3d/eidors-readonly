function show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% Show a 3d view of an object with many slices through it

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: show_3d_slices.m,v 1.2 2007-09-28 01:02:57 aadler Exp $

gca;
hold on

clim = [];
ref_lev= 'use_global';
rimg= calc_slices( img, z_cuts(:)*[inf,inf,1]);
c_img = calc_colours( rimg, clim, 0, ref_lev );
xyz_max= max(img.fwd_model.nodes);
xyz_min= min(img.fwd_model.nodes);

for zi= 1:length(z_cuts)
   surf_slice(rimg(:,:,zi), c_img(:,:,zi), xyz_min, xyz_max, z_cuts(zi), 1);
end

hold off

function surf_slice(rimg, c_img, xyz_min, xyz_max, z_cut, show_surf);
   np= calc_colours('npoints');

   [x,y]=meshgrid( linspace(xyz_min(1),xyz_max(1), np), ...
                   linspace(xyz_min(2),xyz_max(2), np) );

   ff=isnan(rimg);
   c_img(isnan(rimg))= NaN;


   if show_surf
      hh=surf(x,y,c_img*0+z_cut,c_img);
      % WHY WOULD IT BE ANYTHING ELSE  - STUPID MATLAB !!!
      set(hh,'CDataMapping','direct', ...
             'EdgeAlpha',0);
   end

   draw_line_around(c_img, rimg, x,y, [1,0;0,1;0,0], [0;0;z_cut] );


function draw_line_around(c_img, rimg, x,y, M_trans, M_add);
% The MATLAB contour functions are inflexible crap. We
% Need to completely break them to get it to work
% [x_ax;y_ax;z_ax] = M_trans * [x_ax_matlab;y_ax_matlab] + M_add
% For x,z plane at y level 4 we have
% M_trans = [1 0;0 0;0 1]; M_add= [0;4;0];

   c_img(isnan(rimg))= -1e50;
   [jnk,hh]=contour(x,y,c_img, [-1e49,-1e49]);
   Contour_paths = M_trans*[ get(hh,'Xdata'), ...
                             get(hh,'Ydata') ]';
   set(hh,'EdgeColor',[0,0,0], ...
          'Xdata', Contour_paths(1,:) + M_add(1), ...
          'Ydata', Contour_paths(2,:) + M_add(2), ...
          'Zdata', Contour_paths(3,:) + M_add(3));
