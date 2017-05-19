imdl=mk_common_model('d2d4c',16);
clf;
h=show_fem(imdl.fwd_model);
set(h,'FaceColor','none')
set(h,'EdgeColor',[1 1 1]*0.25);
set(h,'LineWidth',1.0);
axis([-0.45 +0.45 0.6 1.2]);
axis off;
hold on;
imdl=mk_common_model('d2d0c',16);
h=show_fem(imdl.fwd_model);
set(h,'FaceColor','none')
set(h,'EdgeColor',[0.1 0.1 1.0]*0.9);
set(h,'EdgeAlpha',0.50);
set(h,'LineWidth',4.0);
hold off;
print('-dpdf','fem.pdf');
