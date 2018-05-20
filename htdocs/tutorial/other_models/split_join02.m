idx  = fmdl1.nodes(:,1)<-0.25;
fmdl1.nodes(idx,1) = -0.25;
fmdl1.nodes(:,1) = fmdl1.nodes(:,1) + 0.25;
subplot(131); show_fem(fmdl1,[0,1,1]); xlim([-1.3,1.3]); axis off

idx  = fmdl2.nodes(:,1)>+0.25;
fmdl2.nodes(idx,1) = +0.25;
fmdl2.nodes(:,1) = fmdl2.nodes(:,1) - 0.25;
subplot(132); show_fem(fmdl2,[0,1,1]); xlim([-1.3,1.3]); axis off

subplot(133); show_fem(fmdl1,[0,1,1]); xlim([-1.3,1.3]);
hold on; hh=show_fem(fmdl2,[0,1,1]); set(hh,'EdgeColor',[0,0,1]);
hold off; axis off

print_convert('split_join02a.png',struct('pagesize',[12,6]));
