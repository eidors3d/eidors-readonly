% GREIT params eval $Id$

for i= 1:length(imgc);
   imgr = imgc{i};
   imgr.calc_colours.npoints = 32;
   params = eval_GREIT_fig_merit(imgr, xyzr);

   fname = sprintf('%s%c', fname_base, i+'a'-1);

   plot(r, params(1,:)); axis([0,0.9,0,2.0]);  ylabel('AR');
   print_convert([fname,'_ar.png'],  '-density 60',0.3);

   plot(r, params(2,:)); axis([0,0.9,-0.15,0.15]);  ylabel('PE');
   print_convert([fname,'_pe.png'],  '-density 60',0.3);

   plot(r, params(3,:)); axis([0,0.9,0,0.4]);  ylabel('RES');
   print_convert([fname,'_res.png'],  '-density 60',0.3);

   plot(r, params(4,:)); axis([0,0.9,0,0.3]);  ylabel('SD');
   print_convert([fname,'_sd.png'],  '-density 60',0.3);

   plot(r, params(5,:)); axis([0,0.9,0,0.6]);  ylabel('RNG');
   print_convert([fname,'_rng.png'],  '-density 60',0.3);
end
