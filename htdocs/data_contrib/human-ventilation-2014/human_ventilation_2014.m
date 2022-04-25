[vv,auxdata]= eidors_readdata('human-ventilation-2014.eit');
vv = real(vv);
FR = 1/(median(diff(auxdata.t_rel))*1e-6);

% Reconstruction Model
   skip5 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};
   fmdl = mk_library_model('adult_male_32el');
   [fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip5{:});
   opt.imgsz = [32 32];
   opt.square_pixels = true;
   opt.noise_figure = 0.5;
   imdl= mk_GREIT_model(mk_image(fmdl,1), 0.20, [], opt);

% Reconstruct
vh = mean(vv(:,790:810),2);
imgr= inv_solve(imdl, vh, vv);
imgs= -calc_slices(imgr); imgs(isnan(imgs))= 0;

figure(1); clf;
axes('position',[0.05,0.5,0.25,0.45]);
imgr.get_img_data.frame_select = 440;
imgr.calc_colours.ref_level = 0;
imgr.calc_colours.backgnd = [1,1,1];
imgr.calc_colours.greylev = 0.1;
show_slices(imgr); axis on; box off;
xposns = 8:3:25;
yposn = 22;
set(gca,'YTick',xposns); ylim([3,30]);
set(gca,'XTick',yposn);
grid

axes('position',[0.33,0.57,0.65,0.30]);
tax = (0:size(imgr.elem_data,2)-1)/FR;
ss = [squeeze(sum(sum(imgs)))/3e2+.25, ...
      squeeze(imgs(xposns,yposn,:))'];
plot(tax,ss,'LineWidth',1);
set(gca,'YTick',[0,0.4]); ylim([-0.1,0.5]);
box off;
cs = cellstr(num2str(xposns(:)));
legend('Sum',cs{:},'Location','NorthWest','Box','off');
