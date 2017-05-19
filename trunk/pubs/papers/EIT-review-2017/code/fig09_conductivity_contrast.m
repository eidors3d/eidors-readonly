[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
space= logspace(-3,3,37);
   vv=[];
i=1;for rat = [.3,1,2,5];
   extra={'ellips',sprintf(['solid ellips = ', ...
          'ellipsoid(0,0,1;0.1,0,0;0,0.1,0;0,0,%f);'],0.1*rat)};
   fmdl= ng_mk_cyl_models(2,[16,1.0],[0.05],extra); 
   fmdl.stimulation = stim; fmdl.meas_select = msel;
   img= mk_image(fmdl,1);
   vh = fwd_solve(img); vh = vh.meas;
   img.elem_data(fmdl.mat_idx{2}) = 2; show_fem(img);

   j=1;for k=space;
      img.elem_data(fmdl.mat_idx{2}) = k;
      vi = fwd_solve(img); vi = vi.meas;
      vv(j,i)= norm(vi-vh);
   j=j+1;end
i=i+1;end
vv = vv ./ (ones(size(vv,1),1)*max(vv,[],1));
co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(gcf,'defaultAxesColorOrder',0.7*co)
semilogx(space,vv);
legend('0.3','1.0','2.0','5.0','Location','NorthEast')
hh=legend('boxoff');
lp=get(hh,'position');
set(hh,'position',lp-[0.35,.10,0,0]);
xlim([min(space), max(space)])
set(gca,'xtick',[1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
set(gcf,'PaperPosition',[2,2,6,3.0]);
box off
print -depsc2 fig09_conductivity_contrast.eps

