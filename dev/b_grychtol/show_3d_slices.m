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

try 
    np = img.calc_colours.npoints;
catch
    np = calc_colours('npoints');
end

%  show_fem(img.fwd_model);
mdl_min = min(img.fwd_model.nodes);
mdl_max = max(img.fwd_model.nodes);
mdl_rng = mdl_max - mdl_min;
np = round(mdl_rng/min(mdl_rng) * np);
[x_cuts, y_cuts, z_cuts] = get_cuts(img,varargin{:});
xvec = linspace(mdl_min(1), mdl_max(1),np(1)+1);
yvec = linspace(mdl_min(2), mdl_max(2),np(2)+1);
zvec = linspace(mdl_min(3), mdl_max(3),np(3)+1);

img.fwd_model.mdl_slice_mapper.x_pts  = xvec(1:end-1) + .5*diff(xvec);
img.fwd_model.mdl_slice_mapper.y_pts  = yvec(1:end-1) + .5*diff(yvec);
for i= 1:length(z_cuts)
    gmdl = mk_grid_model([],xvec,yvec);
    gmdl.nodes(:,3) = z_cuts(i);
    level = [inf inf z_cuts(i)];
    img.fwd_model.mdl_slice_mapper.level = level; 
    pic  = calc_slices(img,level)';
    show_slice(gmdl,pic,img);
    hold on
end

img.fwd_model.mdl_slice_mapper.x_pts  = yvec(1:end-1) + .5*diff(yvec);
img.fwd_model.mdl_slice_mapper.y_pts  = zvec(1:end-1) + .5*diff(zvec);
for i= 1:length(x_cuts)
    gmdl = mk_grid_model([],yvec,zvec);
    gmdl.nodes(:,2:3) = gmdl.nodes;
    gmdl.nodes(:,1) = x_cuts(i);
    level = [x_cuts(i) inf inf];
    img.fwd_model.mdl_slice_mapper.level = level;
    pic  = calc_slices(img,level)';
    show_slice(gmdl,pic,img);
end

img.fwd_model.mdl_slice_mapper.x_pts  = xvec(1:end-1) + .5*diff(xvec);
img.fwd_model.mdl_slice_mapper.y_pts  = zvec(1:end-1) + .5*diff(zvec);
for i= 1:length(y_cuts)
    gmdl = mk_grid_model([],xvec,zvec);
    gmdl.nodes(:,[1 3]) = gmdl.nodes;
    gmdl.nodes(:,2) = y_cuts(i);
    level = [inf y_cuts(i) inf];
    img.fwd_model.mdl_slice_mapper.level = level;
    pic  = calc_slices(img,level)';
    show_slice(gmdl,pic,img);
end
hold off

function show_slice(gmdl,pic,img)
    ff   = find(isnan(pic(:)));
    gmdl.elems([2*ff, 2*ff-1],:)= [];
    gmdl.coarse2fine([2*ff, 2*ff-1],:)= [];
    gmdl.coarse2fine(:,ff)= [];
    gmdl = rmfield(gmdl,'boundary');
    opt.boundary = true;
    gmdl.fwd_model = fix_model(gmdl, opt);
    gim  = mk_image(gmdl,1);
    gim.elem_data(1:2:end) = pic(~isnan(pic));
    gim.elem_data(2:2:end) = pic(~isnan(pic));
    try gim.calc_colours = img.calc_colours; end
    try gim.calc_slices = img.calc_slices; end
    show_fem(gim);

    
function [x_cuts, y_cuts, z_cuts] =  get_cuts(img, varargin)
   mdl_max= max(img.fwd_model.nodes);
   mdl_min= min(img.fwd_model.nodes);
   if nargin==1;
      % Default show 2 z_cuts and 1 x and 1 y cut
       x_cuts= linspace(mdl_min(1), mdl_max(1), 3); x_cuts([1,3])=[];
       y_cuts= linspace(mdl_min(2), mdl_max(2), 3); y_cuts([1,3])=[];
       z_cuts= linspace(mdl_min(3), mdl_max(3), 4); z_cuts([1,4])=[];
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
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl,1); vh= fwd_solve(img);
   load datacom.mat A B;
   img.elem_data(A) = 1.2;
   img.elem_data(B) = 0.8;
   vi = fwd_solve(img);
   imgr = inv_solve(imdl, vh, vi);
   imgr.calc_colours.npoints =64;
   show_3d_slices(img);
   calc_colours('transparency_thresh', 0.25);
   show_3d_slices(imgr,[1.5,2],[],[]);

   subplot(221);  show_3d_slices(imgr,[1,2],0,0.5);
   subplot(222);  show_3d_slices(img ,[1,2],0,0.5);
   
   imgr.calc_colours.transparency_thresh = -1;
   img.calc_colours.transparency_thresh = -1;
   subplot(223);  show_3d_slices(imgr,[1,2],0.3,0.7);
   subplot(224);  
   show_fem(img.fwd_model);
   hold on
   show_3d_slices(img ,[1,2],0.3,0.7);
   axis off
   
   view(10,18);
