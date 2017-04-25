%Repeat with quadratic and cubic finite elements
%Quadratic FEM
fmdl.approx_type    = 'tet10'; %Quadratic
img2 = mk_image(fmdl,1);
img2.fwd_solve.get_all_meas = 1; %Internal voltage
v2 = fwd_solve(img2); 
v2e=v2.meas; v2all=v2.volt;

%Electrode voltages and difference for linear, quadratic and cubic
subplot(211); plot([v1e,v2e,[v2e-v0e]*1]);
legend('1','2','(2-0) x 10','Location','NorthEast')
xlim([1,100]);

print_convert forward_solvers_3d_high_order03a.png
