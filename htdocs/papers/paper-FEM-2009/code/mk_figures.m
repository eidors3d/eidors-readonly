function mk_figures(no)
% Make figures for FEM -errors paper
%
%
% Figure 1: 2D model example

switch no
  case 1; mk_fig_1;
  case 2; mk_fig_2;
  case 3; mk_3d_fig;
  otherwise; error('huh?')
end

function c= contrast;
  c= 0.01;

function imdl= this_mdl;
  imdl = mk_common_model('c2c2',16);
  imdl.RtR_prior = @noser_image_prior;
  imdl.hyperparameter.value = 0.17;

function mk_3d_fig
  set(gcf,'paperposition',[0.25 2.5 6 6*(5/7)]); clf
  imdl= this_mdl;
  maxh=1;
  fmdl= ng_mk_cyl_models([1,1,maxh],[16,0.5],[0.05,0,0.02]); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vh = fwd_solve(img);
  disp([size(fmdl.elems)]);

  extra={'ball','solid ball = cylinder(0.5,0,0;0.5,0,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'}

  [fmdl,mat_idx]= ng_mk_cyl_models( ...
            [1,1,maxh],[16,0.5],[0.05,0,0.02],extra); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  img.elem_data(mat_idx{2}) = 1 + 100*contrast;
  vi = fwd_solve(img);

  axes('position',[0.1,0.5,0.4,0.4]);
  show_fem(img);
  view(-3,60)

  imdl = this_mdl;
  imgr= inv_solve(imdl, vh, vi);
  axes('position',[0.42,0.5,0.4,0.4]);
  set(gca,'YTickLabel',[]);
  show_fem(imgr); axis equal; axis tight

  print -depsc2 ../figures/fig_3d_example.eps

function mk_fig_1
  imdl= this_mdl;
  th= linspace(0,2*pi,50);
  cc= 0.2*cos(th)+0.5; ss= 0.2*sin(th);

  %figure;
  set(gcf,'paperposition',[0.25 2.5 6 6*(5/7)]); clf
  setax= [-0.4,1,-0.5,0.5];

  axes('position',[0.1,0.5,0.4,0.4]);
  [img, vic1, vhc1, n_en] = ng_mdl(imdl, 1, 0)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.5,0.4,0.4]);
  [img, vif1, vhf1, n_en] = ng_mdl(imdl, 0.05, 0)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.1,0.08,0.4,0.4]);
  [img, vic2, vhc2, n_en] = ng_mdl(imdl, 1.00, 1)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.08,0.4,0.4]);
  [img, vif2, vhf2, n_en] = ng_mdl(imdl, 0.05, 1)
  show_fem(img); axis(setax);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  print -depsc2 ../figures/fig_2d_fems.eps

%%% NEXT FIGURE
  set(gcf,'paperposition',[0.25 2.5 6 6]); clf

  axes('position',[0.1,0.5,0.4,0.4]);
  img = inv_solve(imdl,vhc1,vic1); show_fem(img);
  set(gca,'XTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.5,0.4,0.4]);
  img = inv_solve(imdl,vhf1,vif1); show_fem(img);
  set(gca,'XTickLabel',[]);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.1,0.08,0.4,0.4]);
  img = inv_solve(imdl,vhc2,vic1); show_fem(img);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.08,0.4,0.4]);
  img = inv_solve(imdl,vhf2,vif1); show_fem(img);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  imdl.hyperparameter.tgt_data.meas_t1 = vhf1;
  imdl.hyperparameter.tgt_data.meas_t2 = vif1;
  calc_noise_figure(imdl);
  
  print -depsc2 ../figures/fig_2d_reconst.eps



% no obj if rad = 0
function [img, vi, vh, n_en] = ng_mdl(imdl, maxh, rad)
  if rad ==0;
    extra={'',''};
  else
     extra={'ball','solid ball = cylinder(0.5,0,0;0.5,0,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'}
  end

   p = {maxh,rad};
   fmdl = eidors_obj('get-cache', p, 'ng_mk_cyl_models');
   if isempty(fmdl)
      fmdl= ng_mk_cyl_models([0,1,maxh],[16],[0.1,0,0.02],extra); 
      eidors_obj('set-cache', p, 'ng_mk_cyl_models', fmdl);
   end

  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vh = fwd_solve(img);

  pts = interp_mesh(fmdl, 4);
  pts=(pts(:,1,:)-0.5).^2 + (pts(:,2,:)-0).^2;
  img.elem_data = 1 + contrast*mean(pts < 0.2^2,3);
  img.calc_colours.clim = contrast*2;
  vi = fwd_solve(img);
  n_en = [size(fmdl.nodes,1), size(fmdl.elems,1)];


function mk_fig_2
  imdl= this_mdl;
  for i=1:4
     switch i
       case 1; maxh= 1;
       case 2; maxh= 0.05;
       case 3; maxh= 0.02;
       case 4; maxh= 0.01;
       otherwise; error('huh');
     end
     [img, vi1, vh1, n_en] = ng_mdl(imdl, maxh, 0);
     [img, vi2, vh2, n_en] = ng_mdl(imdl, maxh, 1);
     subplot(2,2,i);
     plot([vi1.meas - vh1.meas, vh2.meas - vh1.meas]);
  end



function imdl = assign_mdl(imdl, fmdl);
  imdl.fwd_model.nodes    = fmdl.nodes;
  imdl.fwd_model.elems    = fmdl.elems;
  imdl.fwd_model.gnd_node = fmdl.gnd_node;
  imdl.fwd_model.electrode= fmdl.electrode;
  imdl.fwd_model.boundary = find_boundary(fmdl.elems);
