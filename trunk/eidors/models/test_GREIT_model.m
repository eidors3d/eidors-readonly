xx=[
  -88.5777  -11.4887    4.6893   49.8963  122.7033  150.3033  195.5103 249.7573 ...
  258.8013  279.7393  304.9623  309.2443  322.0923  337.7963  340.6503 348.2633 ...
  357.3043  358.7333  361.5873  364.9183  365.3943  366.3453  366.3453 365.3943 ...
  362.5393  351.5943  343.5053  326.8513  299.2503  288.3073  264.9923 224.0703 ...
  206.4633  162.6833  106.5313   92.2543   57.5153    7.0733   -8.6297 -42.4167 ...
  -90.9547 -105.7057 -134.2577 -178.0367 -193.2647 -222.7687 -265.5957 -278.9197 ...
 -313.1817 -355.5337 -363.6237 -379.3267 -397.8857 -400.7407 -401.6927 -398.8377 ...
 -395.0307 -384.0867 -368.3837 -363.6247 -351.7277 -334.1217 -328.4117 -314.1357 ...
 -291.2947 -282.7297 -267.0257 -236.5707 -221.8187 -196.5977 -159.4807 -147.5837];

yy=[
 -385.8513 -386.8033 -386.3273 -384.8993 -368.7193 -353.9673 -323.0363 -283.5403 ...
 -274.9743 -254.0363 -225.4843 -217.8703 -187.4153 -140.7813 -124.6013  -86.0573 ...
  -38.4703  -29.4273   -9.9173   21.0137   32.4347   53.3727   83.8257   93.3437 ...
  114.7587  149.0237  161.8717  187.5677  222.3037  231.3447  247.5237  267.5087 ...
  271.3177  277.0297  281.3127  279.4097  274.6507  273.2227  276.5547  284.6447 ...
  295.1127  297.4927  301.7757  304.1557  302.2537  297.4947  287.5017  282.2667 ...
  259.9017  225.6387  213.7427  185.6677  141.4127  125.2337   88.5917   34.8187 ...
   17.6897  -22.2803  -73.6723  -85.0923 -117.9263 -163.6083 -176.4573 -205.9613 ...
 -245.9343 -256.4023 -275.4373 -304.9403 -315.4083 -332.0623 -352.0473 -355.3783];

a = [xx; yy]';
a = flipud(a);
a = a - repmat(mean(a),[72 1]); %mk_GREIT_model::stim_targets assumes this
fmdl = ng_mk_extruded_model({300,a,[3,10],25},[16,1.11,150],[1]);
% maxx = max(abs(fmdl.nodes(:,1)));
% maxy = max(abs(fmdl.nodes(:,2)));
% scale = max(maxx,maxy);
% fmdl.nodes = fmdl.nodes/scale;
% fmdl.nodes(:,3) = fmdl.nodes(:,3) - mean(fmdl.nodes(:,3));
fmdl = mdl_normalize(fmdl,1);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;
%%
img = mk_image(fmdl,1);
opt.imgsz = [32 64];
opt.distr = 3; % uniform`
imdl = mk_GREIT_model(img, 0.25, 0.2, opt);

vh=fwd_solve(img);
% figure
% show_fem(img)
select_fcn = inline('(x-20).^2+(y-30).^2+(z-150).^2<10^2','x','y','z');
img.elem_data = 1 + 0.1*elem_select(img.fwd_model, select_fcn);

vi=fwd_solve(img);
figure
show_fem(img);
title('extruded model w/ target')
% print -dpng model.png
figure
imgr= inv_solve(imdl,vh,vi);
title('current GREIT recon')
show_fem(imgr);
axis equal
% print -dpng sln.png

%%
imdl= mk_common_model('c2c2',16);
img=mk_image(imdl); vh=fwd_solve(img); 
img.elem_data(50)=1.1; vi=fwd_solve(img);
figure, show_fem(img), title('circ model');
imgr = inv_solve(imdl,vh,vi);
figure
show_fem(imgr);
title('GN recon')
%%
im_gr = mk_common_gridmdl('GREITc1');
figure
show_fem(inv_solve(im_gr,vh,vi))
title('original GREIT recon');
%%
imdl = mk_common_model('b3cr',[16,1]); %3d, 16 elec, 1 ring
imdl.fwd_model = mdl_normalize(imdl.fwd_model,1);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
imdl.stimulation = stim;
img = mk_image(imdl);
gr_imdl = mk_GREIT_model(img,0.25,0.02);
vh=fwd_solve(img);
img.elem_data(4100)=1.1; 
vi=fwd_solve(img);
figure, show_fem(img), title('cyl model');
figure
show_fem(inv_solve(gr_imdl,vh,vi))
title('current GREIT recon')
axis equal

figure
show_fem(inv_solve(im_gr,vh,vi))
title('original GREIT recon')
axis equal
