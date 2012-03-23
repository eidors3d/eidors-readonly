% ensure dev/m_crabb/forward_problem is on the path

imdl = mk_common_model('f2C0',16);
fmdl = imdl.fwd_model;
img0 = mk_image(fmdl,1);
v0 = fwd_solve(img0); v0=v0.meas;

fmdl.solve = @mc_fwd_solve;
fmdl.system_mat = @mc_calc_system_mat;
fmdl.jacobian   = @mc_calc_jacobian;
fmdl.fem_modify = @mc_fem_modify;
fmdl.mc_type    = 'tri3'; % linear case
img1 = mk_image(fmdl,1);
v1 = fwd_solve(img1); v1=v1.meas;

fmdl.mc_type    = 'tri6'; % linear case
img2 = mk_image(fmdl,1);
v2 = fwd_solve(img2); v2=v2.meas;

fmdl.mc_type    = 'tri20'; % linear case
img3 = mk_image(fmdl,1);
v3 = fwd_solve(img3); v3=v3.meas;

plot([v0,v1,v2,v3,[v1-v0,v2-v0,v3-v0]*100]);
legend('0','1','2','3','1-0','2-0','3-0')
xlim([1,100]);
