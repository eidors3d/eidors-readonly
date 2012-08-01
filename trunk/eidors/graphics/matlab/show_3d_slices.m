function show_3d_slices(img, varargin);
% show_3d_slices(img, z_cuts, x_cuts, y_cuts)
% Show a 3d view of an object with many slices through it
%  z_cuts = planes in z to do a cut
%  x_cuts = planes in x to do a cut
%  y_cuts = planes in y to do a cut
% Default show 2 z_cuts and 1 x and 1 y cut

% (C) 2007-2012 Andy Adler & Bartlomiej Grychtol. 
% License: GPL version 2 or version 3
% $Id$

if isstr(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end
[jnk,ref_lev,max_scale] = scale_for_display( img.elem_data);
try 
    img.calc_colours.ref_level; 
catch
    img.calc_colours.ref_level = ref_lev;
end
try
    img.calc_colours.clim;
catch
    img.calc_colours.clim = max_scale;
end
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


for i= 1:length(z_cuts)
    show_slice(img,[inf inf z_cuts(i)]);
    hold on
end

for i= 1:length(x_cuts)
    show_slice(img,[x_cuts(i) inf inf]);
    hold on
end

for i= 1:length(y_cuts)
    show_slice(img,[inf y_cuts(i) inf]);
    hold on
end
hold off

function show_slice(img, level)
    slc = mdl_slice_mesher(img,level);
    try slc.calc_colours = img.calc_colours; end
    show_fem(slc);

    
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
