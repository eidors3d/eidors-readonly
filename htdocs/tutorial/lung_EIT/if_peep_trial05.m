% Show images $Id: if_peep_trial05.m,v 1.1 2008-07-19 20:58:47 aadler Exp $

for loop = 1:2;
   if loop == 1; img = i_injury; fn= 'a';
   else          img = i_treat;  fn= 'b';
   end

   % image properties
   img.calc_colours.npoints     = 32;
   img.calc_colours.window_range= .5;
   img.calc_colours.ref_level   = 0;
   img.calc_colours.greylev     = 0.01;
   img.calc_colours.backgnd     = [1,1,1];
   img.animate_reconstructions.show_times = 1;

   animate_reconstructions(['if_peep_trial05',fn], img);
end
