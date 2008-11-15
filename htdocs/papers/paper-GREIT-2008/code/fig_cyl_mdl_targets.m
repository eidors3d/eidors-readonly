% Show simulated positions $Id$
load ng_cyl_mdl
load sim_targets

show_fem(fmdl)
crop_model(gca, inline('x-z<-15','x','y','z'))
      view(-90,20)

hold on
[xs,ys,zs]=sphere(10); spclr= [0,.5,.5];
for i=1:1:size(xyzr_pt,2);
   xp=15*xyzr_pt(1,i); yp=15*xyzr_pt(2,i); zp=15*xyzr_pt(3,i); rp=15*xyzr_pt(4,i);
   hh=surf(rp*xs+xp, rp*ys+yp, rp*zs+zp);
   set(hh,'EdgeColor',[0,.4,.4],'FaceColor',[0,.8,.8]);
   hh=text(xp,yp,zp+1,num2str(i)); 
   set(hh,'FontSize',7,'FontWeight','bold','HorizontalAlignment','center');
end
hold off

print -dpng -r100 model_and_target_image.png
