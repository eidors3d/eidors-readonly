% Moving Target
stim=  mk_stim_patterns(32, 1,[0,5],[0,5],{},1);
fmdl= ng_mk_cyl_models([2,1,.1],[32,1.0],[0.05]);
fmdl.nodes(:,3) = fmdl.nodes(:,3) - 1;
fmdl.stimulation = stim;
img = mk_image(fmdl,1);
subplot(211);
show_fem_enhanced(img);

ll= linspace(-0.8,0.8,19);
[vh,vi]= simulate_movement(img, [ll;ll;0*ll;0*ll+0.1]);

imdl = mk_common_model('b2C2',32);
imdl.fwd_model.stimulation = stim;
imdl.hyperparameter.value = .003;
imgr = inv_solve(imdl, vh, vi);
subplot(212);
imgr.show_slices.img_cols = 7;
imgr.calc_colours.clim = .2;
show_slices(imgr);
eidors_colourbar(imgr);


% signal from a range of conductivities
clf
for j=1:2;
  if j==1; str='orthobrick(-0.2,-0.2,0.8;0.2,0.2,1.2);'
  else     str='orthobrick(-0.2,-0.2,0.2;0.2,0.2,1.8);';
  end
stim=  mk_stim_patterns(16, 1,[0,5],[0,5],{},1);
extra={'cube', ['solid cube = ',str]};
fmdl= ng_mk_cyl_models(2,[16,1.0],[0.1],extra); 
fmdl.stimulation = stim;
img= mk_image(fmdl,1);
vh = fwd_solve(img);
subplot(2,2,2*j-1);
img.elem_data(fmdl.mat_idx{2}) = 1.1;
show_fem(img);
cond = 10.^linspace(-5,5,19);
clear vi;
sig = zeros(size(cond));
for i=1:length(cond);
   img.elem_data(fmdl.mat_idx{2}) = cond(i);
   vi(i) = fwd_solve(img);
   sig(i) = mean(abs(vi(i).meas - vh.meas));
   if cond(i)<1; sig(i) = -sig(i); end
end
subplot(2,2,2*j);
plot(log10(cond), sig)
end





return

% change in electrode properties
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
body_geometry.union(1).sphere(1).radius     = 0.1;
body_geometry.union(1).sphere(1).center     = [0 0 0.5];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.08*cos(theta(i)) 1.08*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
    electrode_geometry{i}.keep_material_flag = 1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
img = mk_image(fmdl,1);
stim=  mk_stim_patterns(16, 1,[0,1],[0,1],{},1);
img.fwd_model.stimulation = stim;
vh = fwd_solve(img);

select_fcn = '(x-0.2).^2+(y-0.2).^2<0.1^2';
memb_frac = elem_select( img.fwd_model, select_fcn);
img.elem_data = 1 + memb_frac*(0.1+.05i);
vi = fwd_solve(img);

vi = add_noise(1,vi,vh);


imgs = img;
imgs.elem_data(fmdl.mat_idx{1}) = 1.1;
imgs.elem_data(fmdl.mat_idx{3}) = 1.1;
imgs.fwd_model.electrode = struct([]);
subplot(321);
show_fem(imgs,[1,1]); view(0,80);

subplot(322);
plot([vh.meas,vi.meas,100*(vh.meas-vi.meas)]);


subplot(323);
imdl = mk_common_model('b2C2',16);
imdl.fwd_model.stimulation = stim;
imdl.hyperparameter.value = .003;
imgr = inv_solve(imdl, vh, vi);
imgr.calc_colours.ref_level = 0;
show_fem(imgr,[0,1]);

subplot(325);
%imgr.calc_colours.component = 'imag';
imgr.elem_data = abs( imgr.elem_data );
imgr.calc_colours.ref_level = 0;
show_fem(imgr,[0,1]);


subplot(324);
img.elem_data(fmdl.mat_idx{1}) = 10i;

vi2 = fwd_solve(img);
imgr = inv_solve(imdl, vh, vi2);
imgr.calc_colours.ref_level = 0;
show_fem(imgr,[1,1]);

subplot(326);
imgr.calc_colours.component = 'imag';
imgr.calc_colours.ref_level = 0;
show_fem(imgr,[0,1]);

