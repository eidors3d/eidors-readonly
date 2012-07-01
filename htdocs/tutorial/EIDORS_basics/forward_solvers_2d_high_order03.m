%Repeat with quadratic and cubic finite elements
%Quadratic FEM
fmdl.approx_type    = 'tri6'; %Quadratic
img2 = mk_image(fmdl,1);
img2.fwd_solve.get_all_meas = 1; %Internal voltage
v2 = fwd_solve(img2); 
v2e=v2.meas; v2all=v2.volt;

%Cubic FEM
fmdl.approx_type    = 'tri10'; %Cubic
img3 = mk_image(fmdl,1);
img3.fwd_solve.get_all_meas = 1; %Internal voltage
v3 = fwd_solve(img3); 
v3e=v3.meas; v3all=v3.volt;

%Electrode voltages and difference for linear, quadratic and cubic
figure; plot([v1e,v2e,v3e,[v2e-v0e,v3e-v0e]*10]);
legend('1','2','3','(2-0) x 10','(3-0) x 10')
xlim([1,100]);

print_convert forward_solvers_2d_high_order03a.png
