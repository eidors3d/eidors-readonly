   imdl = mk_common_model('n3r2',[16,2]); fmdl= imdl.fwd_model;
   fmdl1= crop_model(fmdl,inline('x<-0.55','x','y','z'));
   idx  = fmdl1.nodes(:,1)<-0.25;
   fmdl1.nodes(idx,1) = -0.35;
   fmdl1.electrode([19,9])=[];
   fmdl1.nodes(:,1) = fmdl1.nodes(:,1) + 0.35;


   fmdl2 = fmdl1;
   fmdl2.nodes(:,1) = -fmdl2.nodes(:,1);

subplot(121); show_fem(fmdl1); axis off; view(0,65); xlim([-1.3,1.3]);
subplot(122); show_fem(fmdl1); axis off; view(0,65); xlim([-1.3,1.3]);
hold on; hh=show_fem(fmdl2,[0,1,0]); set(hh,'EdgeColor',[0,0,1]);
hold off;

print_convert('split_join04a.png',struct('pagesize',[12,5]));
