function mk_figures(no)
% Make figures for FEM -errors paper
%
%
% Figure 1: 2D model example

switch no
  case 1; mk_fig_1;
  case 1; mk_fig_1;
  otherwise; error('huh?')
end


function mk_fig_1
  imdl = mk_common_model('c2c2',16);
  stim= imdl.fwd_model.stimulation;
  setax= [-0.4,1,-0.7,0.7];
  th= linspace(0,2*pi,50);
  cc= 0.2*cos(th)+0.5; ss= 0.2*sin(th);

  %figure;
  set(gcf,'paperposition',[0.25 2.5 6 6]); clf

  axes('position',[0.1,0.5,0.4,0.4]);
  [img, vi, vh, n_en] = ng_mdl(imdl, 1, 0)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.5,0.4,0.4]);
  [img, vi, vh, n_en] = ng_mdl(imdl, 0.05, 0)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.1,0.08,0.4,0.4]);
  [img, vi, vh, n_en] = ng_mdl(imdl, 1.00, 1)
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.08,0.4,0.4]);
  [img, vi, vh, n_en] = ng_mdl(imdl, 0.05, 1)
  show_fem(img); axis(setax);
  set(gca,'YTickLabel',[]);
  hh=line(cc,ss); set(hh,'Color',[0,0,1],'LineWidth',2); 



% no obj if rad = 0
function [img, vi, vh, n_en] = ng_mdl(imdl, maxh, rad)
  if rad ==0;
    extra={'',''};
  else
     extra={'ball','solid ball = cylinder(0.5,0,0;0.5,0,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'}
  end

  fmdl= ng_mk_cyl_models([0,1,maxh],[16],[0.1,0,0.02],extra); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vh = fwd_solve(img);

  pts = interp_mesh(fmdl, 4);
  pts=(pts(:,1,:)-0.5).^2 + (pts(:,2,:)-0).^2;
  img.elem_data = 1 + 0.1*mean(pts < 0.2^2,3);
  img.calc_colours.clim = 0.2;
  vi = fwd_solve(img);
  n_en = [size(fmdl.nodes,1), size(fmdl.elems,1)];





function imdl = assign_mdl(imdl, fmdl);
  imdl.fwd_model.nodes    = fmdl.nodes;
  imdl.fwd_model.elems    = fmdl.elems;
  imdl.fwd_model.gnd_node = fmdl.gnd_node;
  imdl.fwd_model.electrode= fmdl.electrode;
