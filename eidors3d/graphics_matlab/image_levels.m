function image_levels(img, levels, clim )
% IMAGE_LEVELS(img, levels, clim  ) show slices at levels of an image
% img    = EIDORS image struct
% levels = array of vertical levels
% clim   = colourmap limit (or default if not specified)
%
% $Id: image_levels.m,v 1.4 2005-06-30 10:13:22 aadler Exp $

set(gcf,'NumberTitle','off');
set(gcf,'Name', img.name);
vtx=  img.fwd_model.nodes;
simp= img.fwd_model.elems;
img = img.elem_data;

fc= [];

set_clim= set_colors( img );
if nargin < 3 
    clim= set_clim;
end

for idx= 1:length(levels);
    subplot(2,3,idx);
    lev= levels(idx);
    [fc] = slicer_plot_n(lev,img,vtx,simp);
    view(2);
    grid;
    caxis([-clim,clim]);
    colorbar;
    axis('off');
    title(sprintf('z=%4.2f',lev));
end

function colour_lim = set_colors( sol );

  s= hot(64); s=s(2:60,:);
  s= [flipud(fliplr(s));0,0,0;s]*.8 + .2;
  colormap(s);
  colour_lim= max(abs(sol));
