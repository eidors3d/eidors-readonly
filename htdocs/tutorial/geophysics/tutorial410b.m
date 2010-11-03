% Reconstruct data on Gallery
% $Id$

% homogeneous starting model
background_resistivity= 15.0; % Unit is Ohm.m
background_conductivity= 1./background_resistivity;

gallery_3D_img= eidors_obj('image',gallery_3D_fwd.name);
gallery_3D_img.fwd_model= gallery_3D_fwd;
gallery_3D_img.elem_data= background_conductivity * ...
   ones(size(gallery_3D_img.fwd_model.elems,1),1);

% build the parameter-to-elements mapping
%USE: sparse pilot-point parameterization
sparsity = 1;
%gallery_3D_img= mk_Pilot2DtoFine3D_mapping(gallery_3D_img,sparsity);
gallery_3D_img.fwd_model.coarse2fine = kron(ones(42,1), speye(1024));

gallery_3D_img.rec_model.type = 'fwd_model';
gallery_3D_img.rec_model.name = '2d';
gallery_3D_img.rec_model.elems = gallery_3D_img.fwd_model.misc.model2d.elems;
gallery_3D_img.rec_model.nodes = gallery_3D_img.fwd_model.misc.model2d.nodes;

%disp(['Computing the CC and SS matrices = ' gallery_3D_img.fwd_model.misc.compute_CCandSS]);
%[ref_data,gallery_3D_img]= aa_fwd_solve(gallery_3D_img);
[ref_data]= fwd_solve(gallery_3D_img);
residuals= real_data.meas-ref_data.meas;

%% plot the data
subplot(211);
plot([ref_data.meas,real_data.meas]);
%print -r75 -dpng tutorial410b.png;
