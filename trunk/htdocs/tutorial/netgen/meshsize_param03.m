extra={'ball','solid ball = sphere(7.5,0,5;2);'};
stim = mk_stim_patterns(16,1,'{ad}','{ad}',{},1);

for loop2 = 1:7
   switch loop2;
      case 1; maxh= 2.0;
      case 2; maxh= 1.5;
      case 3; maxh= 1.0;
      case 4; maxh= 0.8;
      case 5; maxh= 0.7;
      case 6; maxh= 0.6;
      case 7; maxh= 0.5;
   end
   fmdl = ng_mk_cyl_models([10,15,maxh],[16,5],[0.5]);
   fmdl.stimulation= stim;
   img = mk_image(fmdl,1);
   vh = fwd_solve(img);

   fmdl = ng_mk_cyl_models([10,15,maxh],[16,5],[0.5],extra);
   fmdl.stimulation= stim;
   img = mk_image(fmdl,1);
   img.elem_data(fmdl.mat_idx{2}) = 1.1;
   vi = fwd_solve(img);

   subplot(121); show_fem(img); view(0,60);

   subplot(122); show_fem(inv_solve( imdl, vh, vi(1))); axis image
   line(xcirc,ycirc,'Color',[0,0.5,0],'LineWidth',2);
   ylabel(sprintf('No. Elems= %d', size(img.fwd_model.elems,1)));

   fname= sprintf('meshsize_param03%c%c.png', 'a', loop2-1+'a');
   print_convert( fname ,'-density 100');
end
