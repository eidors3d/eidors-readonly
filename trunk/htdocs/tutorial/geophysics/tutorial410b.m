% Reconstruct data on Gallery
% $Id: tutorial410b.m,v 1.2 2007-06-15 19:37:45 aadler Exp $

% homogeneous starting model
background_resistivity= 15.0; % Unit is Ohm.m
background_conductivity= 1./background_resistivity;

gallery_3D_img= eidors_obj('image',gallery_3D_fwd.name);
gallery_3D_img.fwd_model= gallery_3D_fwd;
gallery_3D_img.elem_data= background_conductivity * ...
   ones(size(gallery_3D_img.fwd_model.elems,1),1);

% build the parameter-to-elements mapping
%USE: sparse pilot-point parameterization
sparsity = 13;
gallery_3D_img= mk_Pilot2DtoFine3D_mapping(gallery_3D_img,sparsity);

disp(['Computing the CC and SS matrices = ' gallery_3D_img.fwd_model.misc.compute_CCandSS]);
[ref_data,gallery_3D_img]= dg_fwd_solve(gallery_3D_img);
residuals= real_data.meas-ref_data.meas;

%% plot the data
subplot(211);
plot([ref_data.meas,real_data.meas]);
%print -r75 -dpng tutorial410b.png;

