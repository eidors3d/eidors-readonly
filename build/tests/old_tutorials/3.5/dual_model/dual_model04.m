% create fine meshes $Id$

% create base model
mdl_coarse=mk_common_model('a2c0',16);

for model = 1:3
   if model==1
% Model 1: 64 elements
      mdl_str= 'a2c0';
   elseif model==2
% Model 2: 256 elements
      mdl_str= 'b2c0';
   elseif model==3
% Model 3: 576 elements
      mdl_str= 'c2c0';
   end

   mdl_fine= mk_common_model(mdl_str,16);
   mdl_fine.fwd_model.mk_coarse_fine_mapping.n_interp= 150;
   mdl_fine.RtR_prior = @noser_image_prior;
   mdl_fine.hyperparameter.value = 3e-2;

   imdl(model)= mdl_fine;
   imdl(model).fwd_model.coarse2fine = ...
       mk_coarse_fine_mapping( mdl_fine.fwd_model, mdl_coarse.fwd_model);

   subplot(2,3, model)
   show_fem(mdl_fine.fwd_model);
end

print_convert dual_model04a.png
