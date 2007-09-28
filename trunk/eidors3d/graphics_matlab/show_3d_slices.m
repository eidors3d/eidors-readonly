function show_3d_slices(img, varargin);
% show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% Show a 3d view of an object with many slices through it
%  z_cuts = planes in z to do a cut
%  x_cuts = planes in x to do a cut
%  y_cuts = planes in y to do a cut
% Default show 2 z_cuts and 1 x and 1 y cut

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id: show_3d_slices.m,v 1.3 2007-09-28 01:29:31 aadler Exp $

cla;
hold on

[xyz_max, xyz_min, rimg, cimg, ...
          x_cuts, y_cuts, z_cuts] = get_slices(img, varargin);

for zi= 1:length(z_cuts)
   M_trans= [1,0;0,1;0,0];
   M_add= [0,0,z_cuts(zi)];
   surf_slice(rimg(:,:,zi), cimg(:,:,zi), xyz_min, xyz_max, ...
              M_trans, M_add, 1);
end

hold off

function [xyz_max, xyz_min, rimg, cimg, ...
          x_cuts, y_cuts, z_cuts] = get_slices(img, varargin);
   xyz_max= max(img.fwd_model.nodes);
   xyz_min= min(img.fwd_model.nodes);
   if nargin==1;
      % Default show 2 z_cuts and 1 x and 1 y cut
       x_cuts= linspace(xyz_min(3), xyz_max(3), 3); x_cuts([1,3])=[];
       y_cuts= linspace(xyz_min(3), xyz_max(3), 3); y_cuts([1,4])=[];
       z_cuts= linspace(xyz_min(3), xyz_max(3), 4); z_cuts([1,4])=[];
   elseif nargin==2;
       z_cuts= varargin{1}{1};
       x_cuts= [];
       y_cuts= [];
   elseif nargin==3;
       z_cuts= varargin{1}{1};
       y_cuts= varargin{2}{1};
       z_cuts= [];
   elseif nargin==3;
       z_cuts= varargin{1}{1};
       y_cuts= varargin{2}{1};
       z_cuts= varargin{3}{1};
   else 
       error('too many inputs');
   end

   limts= [ z_cuts(:)*[inf,inf,1  ]; ...
            x_cuts(:)*[1  ,inf,inf]; ...
            y_cuts(:)*[inf,1  ,inf]] + 1e-10;

   clim = [];
   ref_lev= 'use_global';
   rimg= calc_slices( img, limts);
   cimg = calc_colours( rimg, clim, 0, ref_lev );

function surf_slice(rimg, cimg, xyz_min, xyz_max, M_trans, M_add, show_surf);
   np= calc_colours('npoints');

   [x,y]=meshgrid( linspace(xyz_min(1),xyz_max(1), np), ...
                   linspace(xyz_min(2),xyz_max(2), np) );

   xyz= reshape([x(:),y(:)]*M_trans', np,np,3);

   ff=isnan(rimg);
   cimg(isnan(rimg))= NaN;


   if show_surf
      hh=surf(xyz(:,:,1)+M_add(1), ...
              xyz(:,:,2)+M_add(2), ...
              xyz(:,:,3)+M_add(3), cimg);
      % WHY WOULD IT BE ANYTHING ELSE  - STUPID MATLAB !!!
      set(hh,'CDataMapping','direct', ...
             'EdgeAlpha',0);
   end

   draw_line_around(cimg, rimg, x,y, M_trans, M_add);


function draw_line_around(cimg, rimg, x,y, M_trans, M_add);
% The MATLAB contour functions are inflexible crap. We
% Need to completely break them to get it to work
% [x_ax;y_ax;z_ax] = M_trans * [x_ax_matlab;y_ax_matlab] + M_add
% For x,z plane at y level 4 we have
% M_trans = [1 0;0 0;0 1]; M_add= [0;4;0];

   cimg(isnan(rimg))= -1e50;
   [jnk,hh]=contour(x,y,cimg, [-1e49,-1e49]);
   Contour_paths = M_trans*[ get(hh,'Xdata'), ...
                             get(hh,'Ydata') ]';
   set(hh,'EdgeColor',[0,0,0], ...
          'Xdata', Contour_paths(1,:) + M_add(1), ...
          'Ydata', Contour_paths(2,:) + M_add(2), ...
          'Zdata', Contour_paths(3,:) + M_add(3));
