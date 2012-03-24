for i=0:2
   img.fwd_model.electrode(1).z_contact=0.1^i; %1,.1,.01
   vh = fwd_solve(img);
   imgc.fwd_model.mdl_slice_mapper.xpts = linspace(-0.25,0.25,200);
   imgc.fwd_model.mdl_slice_mapper.ypts = linspace(0.8,1,100);
   q = show_current(imgc,vh.volt);
   hh=show_fem(imgc);
   set(hh,'EdgeColor',[1,1,1]*.75);
   hold on;
    sy = 1 - (25:5:300)/1000; sx= 0*sy - 0.20;
   hh=streamline(q.xp,q.yp, q.xc, q.yc, sx,sy); set(hh,'Linewidth',2);
    sy = 1 - (25:5:80)/1000; sx= 0*sy - 0.20;
   hh=streamline(q.xp,q.yp,-q.xc,-q.yc,-sx,sy); set(hh,'Linewidth',2);
   title(sprintf('streamlines zc = %5.3f',fmdl.electrode(1).z_contact));
   hold off;
   axis([-.15,.15,0.85,1.02]);

   title(sprintf('current near electrode:  zc = %5.3f',0.1^i));
   print_convert(sprintf('contact_impedance04%c.png','a'+i));
end
