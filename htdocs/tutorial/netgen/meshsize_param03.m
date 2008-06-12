for loop1 = 2
   switch loop1;
      case 1; finelevel='';          loop2max= 7;
      case 2; finelevel='-fine';     loop2max= 7;
      case 3; finelevel='-veryfine'; loop2max= 7;
   end
   for loop2 = 6:loop2max
      switch loop2;
         case 1; maxh= '-maxh=2.0';
         case 2; maxh= '-maxh=1.5';
         case 3; maxh= '-maxh=1.0';
         case 4; maxh= '-maxh=0.8';
         case 5; maxh= '-maxh=0.7';
         case 6; maxh= '-maxh=0.6';
         case 7; maxh= '-maxh=0.5';
      end
      move_the_ball % CALL NETGEN

      subplot(121); show_fem(img); view(90,60);

      subplot(122); show_fem(inv_solve( imdl, vh, vi(1))); axis image
      line(xcirc,ycirc,'Color',[0,0.5,0],'LineWidth',2);
      ylabel(sprintf('No. Elems= %d', size(img.fwd_model.elems,1)));

      fname= sprintf('meshsize_param03%c%c.png', loop1-1+'a', loop2-1+'a');
      print('-dpng','-r100', fname );
   end
end

