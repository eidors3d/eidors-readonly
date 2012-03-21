function [img]=mc_calc_system_mat_point(img)
%Assemble the total stiffness matrix (nothing to be done here)

%Get the forward model
mdl=img.fwd_model;

%Assemble the total system matrix
mdl.solver.At=mdl.solver.Am; 

%Put the forward model back into img
img.fwd_model=mdl;
