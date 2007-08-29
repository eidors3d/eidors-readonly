function image_levels(img, levels, clim )
% IMAGE_LEVELS(img, levels, clim  ) show slices at levels of an image
% img    = EIDORS image struct
% levels = array of vertical levels
%  OR
% levels = [ [z_lev1 ,xpos,ypos], ...
%
% clim   = colourmap limit (or default if not specified)

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: image_levels.m,v 1.15 2007-08-29 09:10:28 aadler Exp $

if exist('OCTAVE_VERSION');
   warning('image_levels does not support octave. Try show_slices');
end

set(gcf,'NumberTitle','off');
set(gcf,'Name', img.name);
fwd_mdl= img.fwd_model;
vtx=  fwd_mdl.nodes;
simp= fwd_mdl.elems;
img_data = img.elem_data;

fc= [];

if size(levels,1)==3
   ll = size(levels,2);
   img_cols = max(levels(2,:));
   img_rows = max(levels(3,:));
else
   ll = length( levels );
   img_cols = ceil( sqrt( ll ));
   img_rows = ceil( ll/ img_cols );
end
   subplot(img_rows,img_cols,1);

% Get geometry Fc
fc = eidors_obj('get-cache', fwd_mdl, 'slicer_plot_fc');
if ~isempty( fc )
    eidors_msg('image_levels: using cached value', 3);
else
   [fc] = slicer_plot_n(levels(1),img_data,vtx,simp);
   eidors_obj('set-cache', fwd_mdl, 'slicer_plot_fc', fc);
   eidors_msg('image_levels: setting cached value', 3);
end

% Set mapped colours
global eidors_colours;
mpc= eidors_colours.mapped_colour;
eidors_colours.mapped_colour = 128;

for idx= 1:length(levels);
    subplot(img_rows,img_cols,idx);
    lev= levels(idx);
    slicer_plot_n(lev,img_data,vtx,simp, fc);
    view(2);
    axis('off');
    axis equal
    title(sprintf('z=%4.2f',lev));
end

% Reset mapped colours
eidors_colours.mapped_colour = mpc;
