function image_levels(img, levels )
% IMAGE_LEVELS(fig_handle, img ) show slices at levels of an image
% fig_handle = handle to a figure
% img = EIDORS3D image struct
% levels = array of vertical levels

set(gcf,'NumberTitle','off');
set(gcf,'Name', img.name);
vtx=  img.fwd_model.nodes;
simp= img.fwd_model.elems;
img = img.elem_data;

fc= [];

for idx= 1:length(levels);
    subplot(2,3,idx);
    lev= levels(idx);
    [fc] = slicer_plot_n(lev,img,vtx,simp);
    view(2);
    grid;
    colorbar;
    axis('off');
    title(sprintf('z=%4.2f',lev));
end
