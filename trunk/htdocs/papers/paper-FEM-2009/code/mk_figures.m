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

  %figure;
  set(gcf,'paperposition',[0.25 2.5 6 6]); clf
  axes('position',[0.1,0.5,0.4,0.4]);
  fmdl= ng_mk_cyl_models(0,[16],[0.1,0,0.02]); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vv.vhc1 = fwd_solve(img);
  pts = interp_mesh(fmdl, 5);
  pts=(pts(:,1,:)-0.5).^2 + (pts(:,2,:)-0).^2;
  img.elem_data = 1 + 0.1*mean(pts < 0.2^2,3);
  img.calc_colours.clim = 0.2;
  vv.vic1 = fwd_solve(img);
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  nv.c1 = [size(fmdl.nodes,1), size(fmdl.elems,1)];

  th= linspace(0,2*pi,50);
  hh=line(0.2*cos(th)+0.5, 0.2*sin(th));
  set(hh,'Color',[0,0,1],'LineWidth',2); 

  axes('position',[0.52,0.5,0.4,0.4]);
  fmdl= ng_mk_cyl_models([0,1,0.05],[16],[0.1,0,0.02]); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vv.vhf1= fwd_solve(img);
  pts = interp_mesh(fmdl, 5);
  pts=(pts(:,1,:)-0.5).^2 + (pts(:,2,:)-0).^2;
  img.elem_data = 1 + 0.1*mean(pts < 0.2^2,3);
  img.calc_colours.clim = 0.2;
  vv.vhf1= fwd_solve(img);
  show_fem(img); axis(setax);
  set(gca,'XTickLabel',[]);
  set(gca,'YTickLabel',[]);
  nv.f1 = [size(fmdl.nodes,1), size(fmdl.elems,1)];

  th= linspace(0,2*pi,50);
  hh=line(0.2*cos(th)+0.5, 0.2*sin(th));
  set(hh,'Color',[0,0,1],'LineWidth',2); 

% extra={'ball','solid ball = cylinder(0.2,0.2,0;0.2,0.2,1;0.2) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.03;'}
  extra={'ball','solid ball = cylinder(0.5,0,0;0.5,0,1;0.2) and orthobrick(-1,-1,0;1,1,0.05);'}

  axes('position',[0.1,0.08,0.4,0.4]);
  fmdl= ng_mk_cyl_models(0,[16],[0.1,0,0.02],extra); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vv.vhc2 = fwd_solve(img);

  ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.5).^2 + (ctr(:,2)-0).^2;
  img.elem_data = 1 + 0.1*(ctr < 0.2^2);
  img.calc_colours.clim = 0.2;
  vv.vic2 = fwd_solve(img);

  show_fem(img); axis(setax);
  th= linspace(0,2*pi,50);
  hh=line(0.2*cos(th)+0.5, 0.2*sin(th));
  set(hh,'Color',[0,0,1],'LineWidth',2); 
  nv.c2 = [size(fmdl.nodes,1), size(fmdl.elems,1)];

  axes('position',[0.52,0.08,0.4,0.4]);
  fmdl= ng_mk_cyl_models([0,1,0.05],[16],[0.1,0,0.02],extra); 
  imdl = assign_mdl(imdl, fmdl);
  img = calc_jacobian_bkgnd(imdl);
  vv.vhf2 = fwd_solve(img);

  ctr = interp_mesh(fmdl); ctr=(ctr(:,1)-0.5).^2 + (ctr(:,2)-0).^2;
  img.elem_data = 1 + 0.1*(ctr < 0.2^2);
  img.calc_colours.clim = 0.2;
  vv.vif2 = fwd_solve(img);

  show_fem(img); axis(setax);
  th= linspace(0,2*pi,50);
  hh=line(0.2*cos(th)+0.5, 0.2*sin(th));
  set(hh,'Color',[0,0,1],'LineWidth',2); 
  set(gca,'YTickLabel',[]);
  nv.f2 = [size(fmdl.nodes,1), size(fmdl.elems,1)];

nv



function imdl = assign_mdl(imdl, fmdl);
  imdl.fwd_model.nodes    = fmdl.nodes;
  imdl.fwd_model.elems    = fmdl.elems;
  imdl.fwd_model.gnd_node = fmdl.gnd_node;
  imdl.fwd_model.electrode= fmdl.electrode;
