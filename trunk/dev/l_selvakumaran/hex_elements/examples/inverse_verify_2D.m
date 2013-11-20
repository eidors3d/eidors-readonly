clc;
clear all;

run EIDORS/startup;

%Import the reference mesh
[ref_mdl, mat_indices]= ls_gmsh_mk_fwd_model('mesh/verify_2D_fine.msh', 'square', 'elec',[], 0.01);
n=size(ref_mdl.elems,1);

%Create the stimulation and measurement patterns
[stim, meas_sel]= mk_plate_stim_patterns(20,1, '{op}', '{op}', [], 1);
ref_mdl.stimulation= stim;
ref_mdl.meas_select= meas_sel;
ref_mdl.normalize_measurements= 1;

%Create the homogeneous image
img_1=mk_image(ref_mdl,1);

%solve the homogeneous problem
img_1.fwd_solve.get_all_meas = 1;
data_1 = fwd_solve(img_1);

%Create the inhomogeneous image
img_2 = img_1;
select_fcn = inline('(x-3).^2+(y-3).^2<1.5^2','x','y','z');
img_2.elem_data = 1 + elem_select(img_2.fwd_model, select_fcn);

%solve the inhomogeneous problem
img_2.fwd_solve.get_all_meas = 1;
data_2 = fwd_solve(img_2);

%Import the working mesh
[fwd_mdl, mat_indices]= ls_gmsh_mk_fwd_model('mesh/verify_2D_coarse.msh', 'square', 'elec',[], 0.01);

%Assign the stimulation and measurement patterns to the fwd_mdl
fwd_mdl.stimulation= stim;
fwd_mdl.meas_select= meas_sel;
fwd_mdl.normalize_measurements= 0;
img_3=mk_image(fwd_mdl,1);

%Create the inverse model
imdl=select_imdl(fwd_mdl);
imdl.fwd_model.jacobian='ls_jacobian_adjoint';
imdl.RtR_prior='ls_prior_laplace';
imdl.solve= 'ls_inv_solve_diff_GN_one_step';
imdl.jacobian_bkgnd.value=2;

% L-curve generation
% hp=[1e-1;1e-2;1e-3;1e-4;1e-5;1e-6;1e-7;1e-8;1e-9;1e-10;1e-11];
% %hp=1e-8;
% x=size(hp);
% y=size(hp);
% for i=1:size(hp)
imdl.hyperparameter.value=1e-8;

%solve the inverse problem
img=inv_solve(imdl,data_1,data_2);
% x(i)=img.residual;
% y(i)=img.normx;
% end
% figure;
% loglog(x,y,'-s');

figure;
h1= subplot(131);
show_fem(img_1,[1,0,0]); 
h2= subplot(132);
show_fem(img_2,[1,0,0]); 
h3= subplot(133);
show_fem(img,[1,0,0]);





