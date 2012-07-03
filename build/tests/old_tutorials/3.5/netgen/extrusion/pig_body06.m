img = imread('pig-thorax.jpg');
colormap(gray(256));
imagesc(.01*([1:371]-438),.01*(-24-[1:371]),img);
set(gca,'YDir','normal');

hh=streamline(q.xp,q.yp, q.xc, q.yc,sx,sy); set(hh,'Linewidth',2);
hh=streamline(q.xp,q.yp,-q.xc,-q.yc,sx,sy); set(hh,'Linewidth',2);

axis equal; axis tight; axis off; print_convert pig_body06a.jpg
