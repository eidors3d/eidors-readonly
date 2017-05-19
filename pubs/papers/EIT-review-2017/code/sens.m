skip_list=[0 0 4];
ne=16;
clf;
for ii=1:length(skip_list)
   ii
   skip_elec=skip_list(ii)
   imdl = mk_common_model('m2C',ne);
   stim=mk_stim_patterns(ne,1,[1 2+skip_elec],[1 2+skip_elec],{},1);
   if ii == 1
      stim = stim_meas_list(stim);
      stim = stim(40,:); % choose an interesting pattern
      stim = stim_meas_list(stim,ne);
   end
   imdl.fwd_model.stimulation = stim;
   img=mk_image(imdl,1);
   J=calc_jacobian(img);
   S=sqrt(sum(J.^2,1)); % sensitivity = column 2-norm of Jacobian
   assert(all(S>=0),'sensitivity is always positive');
   img.elem_data = log10(S(:));
   clf;
   img.calc_colours.cb_shrink_move = [0.4 0.8 0];
   img.calc_colours.ref_level = -4;
   img.calc_colours.clim = 3;
%   img.calc_colours.greylev =-4.5;
   h=show_fem(img,[0 1]);
   set(h,'EdgeColor','none'); axis off;
   print('-dpng',sprintf('sens%d.png',ii));
end
clf; c=eidors_colourbar(img); axis off;
set(c,'Location','SouthOutside','Position',[0.1143 0.1730 0.6750 0.0335]);
tt=get(c,'TickLabels')
L={};
for jj = 1:size(tt,1)
   L{jj} = sprintf('10^{%s}',tt(jj,:));
end
set(c,'TickLabels',L);
print('-dpng','sens-cb.png');
