function show_3d_slices(img, varargin);
% show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% Show a 3d view of an object with many slices through it
%  z_cuts = planes in z to do a cut
%  x_cuts = planes in x to do a cut
%  y_cuts = planes in y to do a cut
% Default show 2 z_cuts and 1 x and 1 y cut

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

cla;
hold on

[xyz_max, xyz_min, rimg, cimg, ...
          x_cuts, y_cuts, z_cuts] = get_slices(img, varargin{:});

for zi= 1:length(z_cuts)
   M_trans= [1,0;0,1;0,0];
   M_add= [0,0,z_cuts(zi)];
   idx= zi;
   surf_slice(rimg(:,:,idx), cimg(:,:,idx), xyz_min, xyz_max, ...
              M_trans, M_add, 1);
end

for xi= 1:length(x_cuts)
   M_trans= [0,0;1,0;0,1];
   M_add= [x_cuts(xi),0,0];
   idx= length(z_cuts) + xi;
   surf_slice(rimg(:,:,idx), cimg(:,:,idx), xyz_min, xyz_max, ...
              M_trans, M_add, 1);
end

for yi= 1:length(y_cuts)
   M_trans= [1,0;0,0;0,1];
   M_add= [0,y_cuts(yi),0];
   idx= length(z_cuts) + length(x_cuts) + yi;
   surf_slice(rimg(:,:,idx), cimg(:,:,idx), xyz_min, xyz_max, ...
              M_trans, M_add, 1);
end

hold off

function [xyz_max, xyz_min, rimg, cimg, ...
          x_cuts, y_cuts, z_cuts] = get_slices(img, varargin);
   xyz_max= max(img.fwd_model.nodes);
   xyz_min= min(img.fwd_model.nodes);
   if nargin==1;
      % Default show 2 z_cuts and 1 x and 1 y cut
       x_cuts= linspace(xyz_min(1), xyz_max(1), 3); x_cuts([1,3])=[];
       y_cuts= linspace(xyz_min(2), xyz_max(2), 3); y_cuts([1,3])=[];
       z_cuts= linspace(xyz_min(3), xyz_max(3), 4); z_cuts([1,4])=[];
   elseif nargin==2;
       z_cuts= varargin{1};
       x_cuts= [];
       y_cuts= [];
   elseif nargin==3;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= [];
   elseif nargin==4;
       z_cuts= varargin{1};
       x_cuts= varargin{2};
       y_cuts= varargin{3};
   else 
       error('too many inputs');
   end

   limts= [ z_cuts(:)*[inf,inf,1  ]; ...
            x_cuts(:)*[1  ,inf,inf]; ...
            y_cuts(:)*[inf,1  ,inf]];

   rimg0= calc_slices( img, limts);
% SURF DOESN'T SHOW THE BLOODY OUTER BOUNDARY, WE NEED TO ADD 4 POINTS
   rimg= NaN*ones(size(rimg0,1)+4,size(rimg0,2)+4,size(limts,1));
   rimg(3:end-2,3:end-2,:)= rimg0;
   cimg = calc_colours( rimg, img);

function xyz= linspace_plus4(lim_min, lim_max, np);
   ooo= ones(length(lim_min),1);
   delta = lim_max(:)-lim_min(:);
   delta = max(delta)*ones(size(delta)); 
   middle=(lim_max(:)+lim_min(:))/2;
   xyz= middle(:)*ones(1,np+4) + ...
        delta/(np-1)*[-(np+3)/2:(np+3)/2]; % for plus4
%       delta/(np-1)*[-(np+1)/2:(np+1)/2]; % for plus2
%       delta/(np-1)*[-(np-1)/2:(np-1)/2]; % for plus0

function surf_slice(rimg, cimg, xyz_min, xyz_max, M_trans, M_add, show_surf);
   np= size(rimg,1)-4;
   lim_min= xyz_min*M_trans;
   lim_max= xyz_max*M_trans;
   xyz= linspace_plus4( lim_min, lim_max, np);
   [x,y]= meshgrid(xyz(1,:),xyz(2,:));
   xyz= reshape( [x(:),y(:)]*M_trans', np+4,np+4,3);

   ff=isnan(rimg);
   bdr= (conv2(double(~ff),ones(3),'same')>0) & ff;
   outbdr = ff & ~bdr;

   ver = eidors_obj('interpreter_version');
   if ver.isoctave ==0 && ver.is64bit
      eidors_msg(['Poor you. You''re using matlab on a 64 bit machine. ' ...
                  'Unfortunately, these versions have serious graphics ' ...
                  'bugs and we can''t show a nice image. Sorry. ' ...
                  'Please bug Mathworks for a fix (good luck!).'],0);
   else
      cimg(outbdr)= NaN;
   end

   if show_surf
      hh=surf(xyz(:,:,1)+M_add(1), ...
              xyz(:,:,2)+M_add(2), ...
              xyz(:,:,3)+M_add(3), flipud(cimg));

      set(hh,'EdgeAlpha',0); % Remove background grid

      % WHY WOULD CDataMapping BE ANYTHING ELSE  - STUPID MATLAB !!!
      % In 64 bit matlab, doing CDataMapping direct on a matrix with
      %   NaN's will crash matlab. Damn.
      set(hh,'CDataMapping','direct');
   end

%  draw_line_around(cimg, rimg, x,y, M_trans, M_add);


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

function do_unit_test
   img = mk_image(mk_common_model('n3r2',16),1);
   show_3d_slices(img,[1,2])
