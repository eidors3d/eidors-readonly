function gallery_gradient(elec_posn,data_tomel,mapping,n_iter,file_prefix)
%
% Driver routine to test the EIDORS functions written for gallery ERT
% Here: tests the GRADIENT algorithm
%
% elec_posn - contains the electrode positions defining the gallery profile
% data_tomel - contains the ERT data written in the "tomel" format
% mapping = 1: uses the element-merging parameterization
%         = 2: uses the sparse pilot-point parameterization
% n_iter: number of iterations
% file_prefix: names of output files begin with file_prefix
%
% Example:
% gallery_gradient(EZG04_Ring1,Data_Ring1_July2004_Wen32_1,2,3,'Test_Results')
%
% Dominique Gibert, June 2007.
%% Definition of model properties
% converts tomel-formatted data to EIDORS structure
real_data= mk_data_tomel(data_tomel,'Mont-Terri data','Wenner protocol');
% defines the conductivity used for the homogeneous starting model
background_resistivity= 15.0; % Unit is Ohm.m
background_conductivity= 1./background_resistivity;

%% create 3D FEM model of the gallery and load homogeneous model
n_rings= 9;
factor= 2;
levels= [-6 -4 -2.5 -1.5 -1 -0.5 -0.25 0 0.25 0.5 1 1.5 2.5 4 6];
gallery_3D_fwd = mk_gallery(elec_posn,data_tomel,n_rings,factor,levels);
figure; show_fem(gallery_3D_fwd); axis square; view(0.,-100.); drawnow;
gallery_3D_img= eidors_obj('image',gallery_3D_fwd.name);
gallery_3D_img.fwd_model= gallery_3D_fwd;
gallery_3D_img.elem_data= ones(size(gallery_3D_img.fwd_model.elems,1),1)*background_conductivity;

%% build the parameter-to-elements mapping
disp('Performing parameters-to-elements mapping');
switch mapping
    case {1} % element-merging parameterization
        disp('parameterization option: elems merging')
        merge_elems= [16 16 16 16 16 16 16 16 ...
            16 16 16 16 16 16 16 16 ...
            32 32 32 32 ...
            32 32 32 32 ...
            128];
        gallery_3D_img= mk_Coarse2DtoFine3D_mapping(gallery_3D_img,merge_elems);
    case {2} % sparse pilot-point parameterization
        disp('Parameterization option: sparse pilot points + interpolation')
        sparsity = 13;
        disp(['Sparsity of pilot points = ' num2str(sparsity)]);
        gallery_3D_img= mk_Pilot2DtoFine3D_mapping(gallery_3D_img,sparsity);
    otherwise
        disp('Wrong parameterization option')
end
disp(['Solving homogeneous 3D model   | background resistivity = ' num2str(background_resistivity)]);
disp(['Computing the CC and SS matrices = ' gallery_3D_img.fwd_model.misc.compute_CCandSS]);
[ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);
residuals= real_data.meas-ref_data.meas;

%% plot the data
disp(['Real data: mean = ' num2str(mean(real_data.meas)) ...
    ' min = ' num2str(min(real_data.meas)) ...
    ' max = ' num2str(max(real_data.meas))]);
disp(['Homogeneous data: mean = ' num2str(mean(ref_data.meas)) ...
    ' min = ' num2str(min(ref_data.meas)) ...
    ' max = ' num2str(max(ref_data.meas))]);
xax= 1:length(ref_data.meas);
figure; axis square; plot(xax,[ref_data.meas,real_data.meas]); drawnow;
figure; axis square; plot(xax,residuals); drawnow;

%% Now attack the inversion
% now: no more necessary to compute time-consuming matrices CC and SS
gallery_3D_img.fwd_model.misc.compute_CCandSS='n';
for k= 1:n_iter
    disp(['Iteration number ',num2str(k)]);
    disp('Begin Jacobian computation');
    jacobian = dg_calc_jacobian(gallery_3D_img);
    disp('End of Jacobian computation');
    [ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);
    residuals= real_data.meas-ref_data.meas;
    disp('Solving the inverse problem');
    svj= svd(jacobian);
    % compute pseudo-inverse using only the largest singular values
    delta_params= pinv(jacobian,svj(1)/20.)*residuals;
    delta_params= delta_params.*gallery_3D_img.params_mapping.perturb;
    gallery_3D_img.params_mapping.params= gallery_3D_img.params_mapping.params + delta_params;
end

%% Solve final model and display results
[ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);
residuals= real_data.meas-ref_data.meas;
save([file_prefix '_results'],'gallery_3D_img','jacobian','real_data','ref_data');

%% plot the results and save figures
disp(['Real data: mean = ' num2str(mean(real_data.meas)) ...
    ' min = ' num2str(min(real_data.meas)) ...
    ' max = ' num2str(max(real_data.meas))]);
disp(['Homogeneous data: mean = ' num2str(mean(ref_data.meas)) ...
    ' min = ' num2str(min(ref_data.meas)) ...
    ' max = ' num2str(max(ref_data.meas))]);
xax= 1:length(ref_data.meas);
fig_data= figure;
    axis square;
    plot(xax,[ref_data.meas,real_data.meas]);
    drawnow;
fig_residuals= figure;
    axis square;
    plot(xax,residuals);
    drawnow;
fig_fem= figure;
    axis square;
    view(30.,80.);
    show_fem(gallery_3D_img);
    drawnow;
fig_slice= figure;
    gallery_3D_img.elem_data= 1./gallery_3D_img.elem_data; % Temporarily convert to resistivity
    show_slices(gallery_3D_img,[inf,inf,0],20,20);
    gallery_3D_img.elem_data= 1./gallery_3D_img.elem_data; % Return to conductivity
    drawnow;
saveas(fig_data,[file_prefix '_data.fig']);
saveas(fig_residuals,[file_prefix '_residuals.fig']);
saveas(fig_fem,[file_prefix '_fem.fig']);
saveas(fig_slice,[file_prefix '_slice.fig']);

end
