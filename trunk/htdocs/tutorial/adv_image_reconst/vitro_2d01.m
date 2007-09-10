% Set eidors default colours
calc_colours('backgnd',[.9,.9,.9]);

% Generate simulation data without noise and standard reconstruction
% Create circular FEM - creates a eidors_mdl type inv_model.
mdl2dim = mk_common_model('b2c');

% Extract measurements from EIT data files and set parameters and
% movement penalty (symbol mu in Soleimani et al. paper).
load(['../../data_contrib/cg_deforming_tank_phantom/oct05/ph1t1/'...
    'avgdata_ph1t1.mat']);

hparameter = 1e-2;
move_vs_conduct = 20;

% Define a eidors_obj Movement model solved by electrode movement
% algorithms.
mdl2dim.hyperparameter.value= hparameter;
mdlM = mdl2dim;
mdlM.fwd_model.conductivity_jacobian = mdlM.fwd_model.jacobian;
mdlM.fwd_model.jacobian = 'aa_e_move_jacobian';
mdlM.RtR_prior = 'aa_e_move_image_prior';
mdlM.aa_e_move_image_prior.parameters = move_vs_conduct;

% Solve inverse problem for mdl2dim and mdlM eidors_obj models.
img2dim = inv_solve(mdl2dim, vvRef, vvAvg1);  % solved no movement algorithms
imgM = inv_solve(mdlM, vvRef, vvAvg1);     % solved with movement algorithms

% ed= img2dim.elem_data;
% lim= .37;
% ed = ed.*(abs(ed)<=lim) + lim*sign(ed).*(abs(ed)>lim);
% img2dim_scl= img2dim;
% img2dim_scl.elem_data= ed;

mysubplot(1,2,1)
show_fem_move(img2dim);
mysubplot(1,2,2)
show_fem_move(imgM, [], 10);

%         % Calculate artefact for each reconstruction
%         load ex_imask.mat;
%         e_space = calc_element_vol(img2dim);
%         artefacts = ~imask.*img2dim.elem_data;
%         amp = sqrt(sum( e_space.*artefacts.^2 ) / sum( e_space ));
%         fprintf('Standard method artefact power = %f \n',amp);
%         tmpM = imgM;
%         tmpM.elem_data = imgM.elem_data(1:256);
%         artefacts = ~imask.*tmpM.elem_data;
%         amp = sqrt(sum( e_space.*artefacts.^2 ) / sum( e_space ));
%         fprintf('Movement method artefact power = %f \n',amp);

