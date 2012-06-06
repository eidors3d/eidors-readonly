% ensure dev/m_crabb/forward_problem is on the path

%Make an inverse model and extract forward model
imdl = mk_common_model('c2C0',16);
fmdl = imdl.fwd_model;

%Make an image and get voltages using eidors default solver
img0 = mk_image(fmdl,1);

%Forward solver using eidors
v0=fwd_solve(img0); v0e=v0.meas;

%Change the solvers
fmdl.solve = @mc_fwd_solve;
fmdl.system_mat = @mc_calc_system_mat;
fmdl.jacobian   = @mc_calc_jacobian;

%Add solver to perform p_refinement
fmdl.fem_modify = @mc_fem_modify;

%Linear, quadratic and cubic
fmdl.mc_type    = 'tri3'; % linear case
img1 = mk_image(fmdl,1);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); v1e=v1.meas;

%Quadratic FEM
fmdl.mc_type    = 'tri6';
img2 = mk_image(fmdl,1);
img2.fwd_solve.get_all_meas = 1; %Internal voltage
v2 = fwd_solve(img2); v2e=v2.meas;

%Cubic FEM
fmdl.mc_type    = 'tri10';
img3 = mk_image(fmdl,1);
v3 = fwd_solve(img3); v3e=v3.meas;

%Electrode voltages
figure; plot([v0e,v1e,v2e,v3e,[v0e-v1e,v2e-v0e,v3e-v0e]*100]);
legend('0','1','2','3','1-0','2-0','3-0')
xlim([1,100]);


%Internal voltage for linear
v1all = v1.volt; 
img1n = rmfield(img1,'elem_data');
img1n.node_data = v1all(1:size(fmdl.nodes,1),1);
figure; show_fem(img1n,[1,0,0]);

%Internal voltage or quadratic
v2all = v2.volt; 
img2n = rmfield(img2,'elem_data');
img2n.node_data = v2all(1:size(fmdl.nodes,1),1);
figure; show_fem(img2n,[1,0,0]);

%Internal voltage difference of linear-quadratic 
img12n=img1n; img12n.node_data=v1all(1:size(fmdl.nodes,1),1)-v2all(1:size(fmdl.nodes,1),1);
figure; show_fem(img12n,[1,0,0]);