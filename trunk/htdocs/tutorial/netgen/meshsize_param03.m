maxh= '-maxh=2.0';
move_the_ball % CALL NETGEN

subplot(121); show_fem(img); view(90,60);

subplot(122); show_fem(inv_solve( imdl, vh, vi(1))); axis image
line(xcirc,ycirc,'Color',[0,0.5,0],'LineWidth',2);

print -dpng -r100 meshsize_param03a.png
