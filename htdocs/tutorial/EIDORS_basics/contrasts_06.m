% $Id$
[stim,msel] = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
space= logspace(-3,3,37);
diam = 0.10;
i=1; vv=[]; for rat = [.3,1,2,5];
   extra={'targ',sprintf(['solid targ = ', ...
          'cylinder(0,0,0;0,0,1;%f) and ', ...
          'plane(0,0,%f;0,0,-1) and plane(0,0,%f;0,0,1);'], ...
           diam,diam*[-1,+1]*rat+1)};
   fmdl= ng_mk_cyl_models([2,1,.1],[16,1.0],[0.05],extra); 
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
clf; subplot(211);

opt.viewpoint = struct('az',-6,'el',13); show_fem_enhanced(img,opt);
print_convert contrasts_06a.jpg

%subplot(212);
semilogx(space,vv,'LineWidth',2);
legend('0.3','1.0','2.0','5.0','Location','SouthEast')
xlim([min(space), max(space)])
set(gca,'xtick',[1e-3,1e-2,1e-1,1,1e1,1e2,1e3]);
print_convert contrasts_06b.png

