%Make common model, and make an image
imdl = mk_common_model('n3r2',16);
fmdl = imdl.fwd_model;
img0 = mk_image(fmdl,1);

%Forward solve using eidors solver (linear)
v0 = fwd_solve(img0); v0e=v0.meas;

%Solvers (including mc_fem_modify and element type)
fmdl.solve = @mc_fwd_solve;
fmdl.system_mat = @mc_calc_system_mat;
fmdl.jacobian   = @mc_calc_jacobian;
fmdl.fem_modify = @mc_fem_modify;

%Forward solve using mc solvers (linear, quadratic)

%Linear FEM
fmdl.mc_type    = 'tet4';
img1 = mk_image(fmdl,1);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); v1e=v1.meas;

%Quadratic FEM
fmdl.mc_type    = 'tet10';
img2 = mk_image(fmdl,1);
img2.fwd_solve.get_all_meas = 1; %Internal voltage
v2 = fwd_solve(img2); v2e=v2.meas;

%Electrode voltages
figure; plot([v0e,v1e,v2e,[v0e-v1e,v2e-v0e]*2]);
legend('0','1','2','1-0','2-0')
xlim([1,100]);


%Internal voltage for linear
v1all = v1.volt; 
img1n = rmfield(img1,'elem_data');
img1n.node_data = v1all(1:size(fmdl.nodes,1),1);
figure; show_slices(img1n,[inf,inf,2.5]);

%Internal voltage or quadratic
v2all = v2.volt; 
img2n = rmfield(img2,'elem_data');
img2n.node_data = v2all(1:size(fmdl.nodes,1),1);
figure; show_slices(img2n,[inf,inf,2.5]);

%Internal voltage difference of linear-quadratic 
img12n=img1n; img12n.node_data=v1all(1:size(fmdl.nodes,1),1)-v2all(1:size(fmdl.nodes,1),1);
figure; show_slices(img12n,[inf,inf,2.5]);