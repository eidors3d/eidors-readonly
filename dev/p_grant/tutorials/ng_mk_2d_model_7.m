% circle
th = linspace(0,2*pi,26)'; th(end)=[];
circle = 0.1 + 0.2*[cos(th), sin(th)];

% electrode positions
ep = [0 0.4; 0.5 0.4; -0.2 -0.75];

% electrode shapes
es = [.1 15; % small electrode
      0  10; % point electrode
      .2 15];% large electrode

% model with a hole
mdl1 = ng_mk_2d_model({xy,circle,0.05},ep,es);
show_fem(mdl1)
hold on

% boundary of the hole
bb = unique(find_boundary(mdl1));
idx = mdl1.nodes(bb,1) < 0.31 & ...
      mdl1.nodes(bb,1) > -0.11& ...
      mdl1.nodes(bb,2) < 0.31 & ...
      mdl1.nodes(bb,2) > -0.11;
bb = bb(idx);

% order the points counter-clockwise
bp = order_loop(mdl1.nodes(bb,:),-1);

% change netgen options to discourage additional nodes on contours
ng_write_opt('meshoptions.fineness',1,'options.meshsize',0.05);

% model to fill the hole
mdl2 = ng_mk_2d_model(bp);
hh = show_fem(mdl2);
set(hh, 'EdgeColor','red');
hold off

% merge the two models
mdl = merge_meshes(mdl1,mdl2);

% make an image
img = mk_image(mdl,1);
img.elem_data(mdl.mat_idx{2}) = 0.3;
show_fem(img,[0 1 0])

% remove the ng.opt file
delete('ng.opt');