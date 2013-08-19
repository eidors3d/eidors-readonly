ci = [5 .5 .05]; 
for i=1:3
   img.fwd_model.electrode(1).z_contact=ci(i);
   vh = fwd_solve(img);
   imgc.fwd_model.mdl_slice_mapper.xpts = linspace(-0.25,0.25,200);
   imgc.fwd_model.mdl_slice_mapper.ypts = linspace(0.8,1,100);
   q = show_current(imgc,vh.volt);
   hh=show_fem(imgc);
   set(hh,'EdgeColor',[1,1,1]*.75);
   hold on;

   sy = linspace(.98,.8 ,20); sx= 0*sy - 0.15;
   hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
   hh=streamline(q.xp,q.yp,-q.xc,-q.yc,-sx,sy); set(hh,'Linewidth',2);

   title(sprintf('streamlines zc = %5.3f',fmdl.electrode(1).z_contact));
   hold off;
   axis([-.15,.15,0.85,1.02]);

   title(sprintf('current near electrode:  zc = %5.3f',ci(i)));
   print_convert(sprintf('contact_impedance04%c.png','a'+i-1));
end
