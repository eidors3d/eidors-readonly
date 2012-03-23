% ensure dev/m_crabb/forward_problem is on the path

imdl = mk_common_model('n3r2',16);
fmdl = imdl.fwd_model;
img0 = mk_image(fmdl,1);
v0 = fwd_solve(img0); v0=v0.meas;

fmdl.solve = @mc_fwd_solve;
fmdl.system_mat = @mc_calc_system_mat;
fmdl.jacobian   = @mc_calc_jacobian;
fmdl.fem_modify = @mc_fem_modify;
fmdl.mc_type    = 'tet4'; % linear case
img1 = mk_image(fmdl,1);
v1 = fwd_solve(img1); v1=v1.meas;

fmdl.mc_type    = 'tet10'; % linear case
img2 = mk_image(fmdl,1);
img2.fwd_solve.get_all_meas = 1;
v2 = fwd_solve(img2); v2all = v2.volt; v2=v2.meas;
img2n = rmfield(img2,'elem_data');
vtxidx = 1:size(fmdl.nodes,1);
img2n.node_data = v2all(vtxidx,1);


plot([v0,v1,v2,[v1-v0,v2-v0]*1]);
legend('0','1','2','1-0','2-0');
xlim([1,100]);
