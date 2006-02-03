function fig_figures
% Generate fingers images for EIDORS paper
% $Id: fig_fingers.m,v 1.2 2006-02-03 22:19:31 aadler Exp $

calc_colours('mapped_colour', 128); % so matlab can print to eps properly
calc_colours('greylev',0.1);
calc_colours('sat_adj',0.98);
calc_colours('backgnd',[1,1,1]);

imdl=mk_common_model('n3r2');
im=eidors_obj('image','3 fingers','fwd_model',imdl.fwd_model);
rr= ones(size(imdl.fwd_model.elems,1),1);
im.elem_data=rr;
vh= fwd_solve(im);

fingerelems=  [ ...
 232 233 234 238 239 240 253 254 255 ...
 508 509 510 514 515 516 529 530 531];
rr(fingerelems)=1.5;
im.elem_data=rr;
vi= fwd_solve(im);
calc_colours('ref_level',1);
% show_fem(im)

imdl.hyperparameter.value= 1e-4;
irec= inv_solve(imdl,vh,vi);
calc_colours('ref_level',0);
% show_fem(irec)
% show_slices(irec,4)
levels= [.1,.83,1.1,1.72,2.1,2.63];
ll= length(levels);
for i=1:ll
   subplot(2,ll,i);
   slicer_plot(irec, levels(i));
end  

function slicer_plot(img,level)
   fwd_mdl= img.fwd_model;
   vtx=  fwd_mdl.nodes;
   simp= fwd_mdl.elems;
   img_data = img.elem_data;

   fc = eidors_obj('get-cache', fwd_mdl, 'slicer_plot_fc');
   if ~isempty( fc )
       eidors_msg('image_levels: using cached value', 3);
   else
      [fc] = slicer_plot_n(level,img_data,vtx,simp);
      eidors_obj('set-cache', fwd_mdl, 'slicer_plot_fc', fc);
      eidors_msg('image_levels: setting cached value', 3);
   end

   slicer_plot_n(level,img_data,vtx,simp, fc);
   clim = max(abs(img_data));
   caxis([-clim,clim]);
   axis('square'); axis('off');
