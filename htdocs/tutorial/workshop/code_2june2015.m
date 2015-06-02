
imb=  mk_common_model('a2c2',16);
img= mk_image(imb,1);
show_fem(img,[0,0,1]);
vv = fwd_solve(img); vv = vv.meas;
p1 = [54,55];
p2 = [41,40];
sp = logspace(-.7,.7,20);
for i = length(sp):-1:1
   for j = length(sp):-1:1
      img.elem_data(p1) = sp(i);
      img.elem_data(p2) = sp(j);
      v2 = fwd_solve(img); 
      M(i,j) = norm(vv-v2.meas);
   end
end
colormap(gray(256));
imagesc(M);

return

imb=  mk_common_model('c2c',16);
e= size(imb.fwd_model.elems,1);
bkgnd= 1;

% Solve Homogeneous model
img= mk_image(imb.fwd_model, bkgnd);
show_fem(img);

%Determine homogeneous voltages
vh= fwd_solve( img );

% Add Two triangular elements and remake a model
img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])=bkgnd * 2;
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=bkgnd * 2;
img_i = img;
 show_fem(img_i);

%Solve inhomogeneous model
vi= fwd_solve( img_i );


%Add some noise -18dB SNR
vi_n= vi; 
nampl= std(vi.meas - vh.meas)*10^(-18/20);
vi_n.meas = vi.meas + nampl *randn(size(vi.meas));

%Plot voltages with and with
 plot(vi_n.meas,'r'); hold on; plot(vi.meas,'b'); plot(100*(vi_n.meas-vi.meas),'g');

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

%AVOID INVERSE CRIME
% Create another coarser Inverse Model
%Supply it's forward model to inverse model created
imb2=  mk_common_model('b2c',16);
inv2d.fwd_model= imb2.fwd_model;



%% RECONSTRUCTION

%% Gauss-Newton solvers - 
inv2d.solve=       @inv_solve_diff_GN_one_step;


% Tikhonov prior
inv2d.hyperparameter.value = .003;
inv2d.RtR_prior=   @prior_laplace;
imgr1= inv_solve( inv2d, vh, vi);
imgn1= inv_solve( inv2d, vh, vi_n);

subplot(331); show_slices(img_i); 
title LAPLACE
subplot(334); show_slices(imgr1);
subplot(337); show_slices(imgn1);

%NOSER prior
inv2d.hyperparameter.value = .01;
inv2d.RtR_prior=   @prior_noser;
imgr2= inv_solve( inv2d, vh, vi);
imgn2= inv_solve( inv2d, vh, vi_n);

subplot(332); show_slices(img_i); 
title NOSER
subplot(335); show_slices(imgr2); 
subplot(338); show_slices(imgn2);


inv2d.hyperparameter.value = 1e-0;
inv2d.solve=       @inv_solve_diff_pdipm;
inv2d.R_prior=     @prior_TV;
inv2d.parameters.max_iterations= 10;
inv2d.inv_solve_diff_pdipm.norm_data=1;
inv2d.inv_solve_diff_pdipm.norm_prior=1;
imgr3= inv_solve( inv2d, vh, vi);
imgn3= inv_solve( inv2d, vh, vi_n);
subplot(333); show_slices(img_i); 
title PDIPM
subplot(336); show_slices(imgr3); 
subplot(339); show_slices(imgn3);

fmdl = ng_mk_cyl_models([2,1,0.1],[16,1],[0.05]); 
fmdl.stimulation = imgr3.fwd_model.stimulation;
opt.distr = 0; % best for cylinders
i_grc = mk_GREIT_model(fmdl,.2,.01,opt);

imgr3= inv_solve( i_grc, vh, vi);
imgn3= inv_solve( i_grc, vh, vi_n);
subplot(333); show_slices(img_i); 
title GREIT
subplot(336); show_slices(imgr3); 
subplot(339); show_slices(imgn3);

return

e= size(imb.fwd_model.elems,1);
bkgnd= 1;

% Solve Homogeneous model
img= mk_image(imb.fwd_model, bkgnd);
show_fem(img);

%Determine homogeneous voltages
vh= fwd_solve( img );

% Add Two triangular elements and remake a model
img.elem_data([25,37,49:50,65:66,81:83,101:103,121:124])=bkgnd * (1+1i);
img.elem_data([95,98:100,79,80,76,63,64,60,48,45,36,33,22])=bkgnd * (1-1i);
img_i = img;
img_i.calc_colours.component = 'imag';
 show_fem(img_i);

%Solve inhomogeneous model
vi= fwd_solve( img_i );

%Plot voltges and differences
subplot(211);
 plot(vh.meas,'r'); hold on; plot(vi.meas,'g'); plot((vi.meas-vh.meas),'b');

%Plot voltages and relative differences
subplot(212);
 plot(vh.meas,'r'); hold on; plot(vi.meas,'g'); plot((vi.meas-vh.meas)./vh.meas,'b');


%Add some noise -18dB SNR
vi_n= vi; 
nampl= std(vi.meas - vh.meas)*10^(-18/20);
vi_n.meas = vi.meas + nampl *randn(size(vi.meas));

%Plot voltages with and with
 plot(vi_n.meas,'r'); hold on; plot(vi.meas,'b'); plot(100*(vi_n.meas-vi.meas),'g');

% Create Inverse Model
inv2d= eidors_obj('inv_model', 'EIT inverse');
inv2d.reconst_type= 'difference';
inv2d.jacobian_bkgnd.value= 1;

%AVOID INVERSE CRIME
% Create another coarser Inverse Model
%Supply it's forward model to inverse model created
imb2=  mk_common_model('b2c',16);
inv2d.fwd_model= imb2.fwd_model;

%Plot the two model
 subplot(121); show_fem(imb.fwd_model); 
subplot(122); show_fem(imb2.fwd_model);


%% RECONSTRUCTION

%% Gauss-Newton solvers - 
inv2d.solve=       @inv_solve_diff_GN_one_step;

% Tikhonov prior
inv2d.hyperparameter.value = .003;
inv2d.RtR_prior=   @prior_tikhonov;
imgr1= inv_solve( inv2d, vh, vi);
imgn1= inv_solve( inv2d, vh, vi_n);
imgr1.calc_colours.component = 'imag';
imgn1.calc_colours.component = 'imag';

subplot(331); show_slices(img_i); 
subplot(334); show_slices(imgr1);
subplot(337); show_slices(imgn1);

% Tikhonov prior
inv2d.hyperparameter.value = .003;
inv2d.RtR_prior=   @prior_laplace;
imgr1= inv_solve( inv2d, vh, vi);
imgn1= inv_solve( inv2d, vh, vi_n);
imgr1.calc_colours.component = 'imag';
imgn1.calc_colours.component = 'imag';

subplot(331); show_slices(img_i); 
title LAPLACE
subplot(334); show_slices(imgr1);
subplot(337); show_slices(imgn1);

%NOSER prior
inv2d.hyperparameter.value = .01;
inv2d.RtR_prior=   @prior_noser;
imgr2= inv_solve( inv2d, vh, vi);
imgn2= inv_solve( inv2d, vh, vi_n);
imgr2.calc_colours.component = 'imag';
imgn2.calc_colours.component = 'imag';

subplot(332); show_slices(img_i); 
title NOSER
subplot(335); show_slices(imgr2); 
subplot(338); show_slices(imgn2);


inv2d.hyperparameter.value = 1e-0;
inv2d.solve=       @inv_solve_diff_pdipm;
inv2d.R_prior=     @prior_TV;
inv2d.parameters.max_iterations= 10;
inv2d.inv_solve_diff_pdipm.norm_data=1;
inv2d.inv_solve_diff_pdipm.norm_prior=1;
imgr3= inv_solve( inv2d, vh, vi);
imgn3= inv_solve( inv2d, vh, vi_n);
imgr3.calc_colours.component = 'imag';
imgn3.calc_colours.component = 'imag';
subplot(333); show_slices(img_i); 
title PDIPM
subplot(336); show_slices(imgr3); 
subplot(339); show_slices(imgn3);



img_i.calc_colours.component = 'imag';

return












return

imdl = mk_common_model('d2d2c',16);
img = mk_image(imdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);
J0= calc_jacobian( img);
sel = elem_select( img.fwd_model,...
  ['(2*(x-0.4).^2 + (y-0.2).^2 < 0.5^2)' ...
  '|(2*(x+0.4).^2 + (y-0.2).^2 < 0.5^2)']);
img.elem_data = img.elem_data - ...
   sel*0.1; 
J1= calc_jacobian( img);

subplot(221); show_fem(img);
sens = img;
vol = get_elem_volume(img.fwd_model);
sens.calc_colours.clim = .5;
subplot(222)
sens.elem_data = J0(1,:)'./vol;
show_fem(sens);
subplot(223)
sens.elem_data = J1(1,:)'./vol;
show_fem(sens);
subplot(224)
sens.calc_colours.clim = .05;
sens.elem_data=(J0(1,:)-J1(1,:))'./vol;
show_fem(sens);
return

fmdl = ng_mk_cyl_models( ...
  [1,1,.10], [16,.5],[0.05]);
img = mk_image(fmdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);
J = calc_jacobian( img);
sens = img;
vol = get_elem_volume(fmdl);
sens.elem_data = J(1,:)'./vol;
show_3d_slices(sens,[0.5],[0],[0])
sens.calc_colours.clim = 0.2;
show_slices(sens,[inf, inf, 0.5]);

is = calc_slices(sens,[inf, inf, 0.5]);
disp([max(abs(is(:))), min(abs(is(:)))]);
disp([max(abs(is(:)))/min(abs(is(:)))]);

return
imdl = mk_common_model('b2c2',16);
img = mk_image(imdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);
J = calc_jacobian( img);
sens = img;
sens.elem_data = J(1,:)';
show_fem(sens)

return
imdl = mk_common_model('b2c2',16);
img = mk_image(imdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);
vh = fwd_solve(img);
sens = img;

Delta_c = 1e-7;
for i=1:length(img.elem_data);
  imgt = img; imgt.elem_data(i)= ...
     imgt.elem_data(i) + Delta_c;
  vi = fwd_solve(imgt);
  sens.elem_data(i) = ...
    (vi.meas - vh.meas)/Delta_c; 
end
show_fem(sens);

return
switch 1;
  case 1;
   imdl = mk_common_model('b2c2',16);
   fmdl = imdl.fwd_model;
  case 2;
   fmdl = ng_mk_cyl_models( ...
     [1,1,.10], [16,.5],[0.05,0,.005]);
end

img = mk_image(fmdl,1);
show_fem(img,[0,1]);
stim= mk_stim_patterns( ...
    16,1,[0,1],[0,1],{},16);
img.fwd_model.stimulation = stim(1:2);
vh = fwd_solve(img);
plot(vh.meas);

return
load neonate_2006
clf; plot(neonate_p04p_1016)

return
imdl = mk_common_model('b2c2',16);
img = mk_image(imdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);

vh = fwd_solve(img);
show_fem(img);

idx = linspace(-20,20,50);
for i = 1:length(idx);
   img.elem_data([42,55])=1 + 1i*idx(i);
   vi = fwd_solve(img);
   d(i) = vi.meas - vh.meas;
end
subplot(211)
plot(idx,real(d),'LineWidth',3);
subplot(212)
plot(idx,imag(d),'LineWidth',3);





return
imdl = mk_common_model('b2c2',16);
img = mk_image(imdl,1);
show_fem(img,[0,1]);
img.fwd_model.stimulation = ...
   stim_meas_list([12,14,4,6],16);

vh = fwd_solve(img);
show_fem(img);

idx = logspace(-3,3,20);
for i = 1:length(idx);
   img.elem_data([42,55])= idx(i);
   vi = fwd_solve(img);
   d(i) = vi.meas - vh.meas;
end
semilogx(idx,d,'LineWidth',3);
