% Show simulated positions
load ng_mdl_16x1_fine;
show_fem(ng_mdl_16x1_fine)
crop_model(gca, inline('x-z<-15','x','y','z'))
      view(-90,70)
%     view(-90,20)

hold on
[xs,ys,zs]=sphere(10); spclr= [0,1,1];
for i=1:5:size(xyzr_pt,2);
   xp=xyzr_pt(1,i);
   yp=xyzr_pt(2,i);
   zp=xyzr_pt(3,i);
   rp=xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',spclr,'FaceColor',spclr);
end
hold off

% print -dpng -r100 simulation_test02a.png

