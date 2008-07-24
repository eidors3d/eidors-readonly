function [im0,irec]= fig_figures
% Generate fingers images for EIDORS paper
% $Id: fig_fingers.m,v 1.4 2006-02-07 03:54:39 aadler Exp $

calc_colours('mapped_colour', 128); % so matlab can print to eps properly
calc_colours('greylev',0.1);
calc_colours('sat_adj',0.95);
calc_colours('backgnd',[1,1,1]);
calc_colours('ref_level',1);

imdl=mk_common_model('n3r2');
im0=eidors_obj('image','3 fingers','fwd_model',imdl.fwd_model);
rr= ones(size(imdl.fwd_model.elems,1),1);
im0.elem_data=rr;
vh= fwd_solve(im0);

fingerelems=  [ ...
 232 233 234 238 239 240 253 254 255 ...
 508 509 510 514 515 516 529 530 531];
rr(fingerelems)=1.5;
im0.elem_data=rr;
vi= fwd_solve(im0);
% show_fem(im0)

imdl.hyperparameter.value= 1e-4;
irec= inv_solve(imdl,vi,vh);
irec.elem_data= irec.elem_data + calc_colours('ref_level');

%slicer_plots(im0 ,'simulated_inhomogeneities.eps');
%slicer_plots(irec,'reconstructed_conductivity.eps');
   clf;
 slicer_plots(im0 ,'', .35);
 slicer_plots(irec,'fig_fingers.eps', 0);

function slicer_plots(img,fname,tb)
   levels= [.1,.83,1.1,1.72,2.1,2.63];
   ll= length(levels);
   dx= 0.9/ll;
   for i=1:ll
%     subplot(1,ll,i);
      axes('position', [.05+dx*(i-1),.05+tb,dx*1.1,.4 ]);
      slicer_plot(img, levels(i));
      % kill colorbar. We must create then kill the colorbar,
      %   otherwise the figures with and without it look different.
      if i<ll
         hh= colorbar; axis(hh,'off');
      end
   end  

   if nargin>=2 && ~isempty(fname)
      set(gcf,'PaperUnits','inches');
      set(gcf,'PaperPosition',[.25,2.5,8,4]); 
      print(gcf,'-depsc2',fname);
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
