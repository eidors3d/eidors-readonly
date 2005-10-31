function image_levels(img, levels, clim )
% IMAGE_LEVELS(img, levels, clim  ) show slices at levels of an image
% img    = EIDORS image struct
% levels = array of vertical levels
% clim   = colourmap limit (or default if not specified)

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: image_levels.m,v 1.9 2005-10-31 01:21:12 aadler Exp $

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

set_clim= set_colors( img_data );
if nargin < 3 
    clim= set_clim;
end

ll = length( levels );
img_cols = ceil( sqrt( ll ));
img_rows = ceil( ll/ img_cols );
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

for idx= 1:length(levels);
    subplot(img_rows,img_cols,idx);
    lev= levels(idx);
    slicer_plot_n(lev,img_data,vtx,simp, fc);
    view(2);
    grid;
    caxis([-clim,clim]);
    colorbar;
    axis('off');
    title(sprintf('z=%4.2f',lev));
end

function colour_lim = set_colors( sol );
   global eidors_colours;
   mpc= eidors_colours.mapped_colour;
   eidors_colours.mapped_colour = 0;
   colormap( squeeze( calc_colours( linspace(-1,1,128) ) ));
   colour_lim= max(abs(sol));
   eidors_colours.mapped_colour = mpc;
