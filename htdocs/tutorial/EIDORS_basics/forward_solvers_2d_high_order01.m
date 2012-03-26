%Make an inverse model and extract forward model
imdl = mk_common_model('c2C0',16);
fmdl = imdl.fwd_model;


%Default EIDORS solver
%Make image of unit conductivity
img0 = mk_image(fmdl,1);
img0.fwd_solve.get_all_meas = 1; %Internal voltage
v0=fwd_solve(img0);
v0e=v0.meas; v0all=v0.volt; 


%High-order EIDORS solver
%Change default eidors solvers
fmdl.solve = @mc_fwd_solve;
fmdl.system_mat = @mc_calc_system_mat;
%Add extra solver for polynomial refinement and add element type
fmdl.fem_modify = @mc_fem_modify; fmdl.mc_type    = 'tri3'; % linear
%Make an image and get voltages using high order solver
img1 = mk_image(fmdl,1);
img1.fwd_solve.get_all_meas = 1; %Internal voltage
v1 = fwd_solve(img1); 
v1e=v1.meas; v1all=v1.volt;


%Plot electrode voltages and difference
figure; plot([v0e,v1e,[v0e-v1e]*100]);
legend('0','1','(1-0) x 100'); xlim([1,100]);

print_convert forward_solvers_2d_high_order01a.png


