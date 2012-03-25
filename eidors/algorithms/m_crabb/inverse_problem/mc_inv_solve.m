function img= mc_inv_solve( inv_model, data1, data2)
% MC_INV_SOLVE inverse solver 
% img= mc_INV_solve( inv_model, data1, data2)
% img        => output image (or vector of images)
% inv_model  => inverse model struct
% data1      => differential data at earlier time
% data2      => differential data at later time
%
% both data1 and data2 may be matrices (MxT) each of
%  M measurements at T times
% if either data1 or data2 is a vector, then it is expanded
%  to be the same size matrix

%Calculate the voltage difference data
diff_volt = calc_difference_data( data1, data2, inv_model.fwd_model);

%Get the change in conductivity from generalised tikhonov regularisation
sol = tikhonov_one_step_inv(inv_model) * diff_volt;

%Create a data structure to return (difference conductivity!!!!)
img.name= 'solved by mc_GN_solve';
img.elem_data = sol;
img.fwd_model= inv_model.fwd_model;
img.type='image';

function tikhonov_inv = tikhonov_one_step_inv(inv_model)
% The one_step reconstruction matrix is cached
   tikhonov_inv = eidors_obj('get-cache', inv_model, 'mc_2011_one_step_inv');
   if ~isempty(tikhonov_inv)
       eidors_msg('mc_inv_solve: using cached value', 3);
   else
       %Calculate the background for the Jacobian
       %For example : imdl.jacobian_bkgnd.value=img_prior.elem_data
       img_bkgnd= calc_jacobian_bkgnd( inv_model );
       %Calculate the Jacobian
       %To make sense for p-refine, need fwd_model.mc_type!!!
       J = calc_jacobian( inv_model.fwd_model, img_bkgnd);
       
       %Calculate the approximation to pd operator for image prior matrix
       %For example : imdl.RtR_prior=@tikhonob_image_prior
       RtR = calc_RtR_prior( inv_model );
       %Calculate the hyperparameter
       %For example : imdl.hyperparameter.value=1e-3(default)
       hp= calc_hyperparameter( inv_model );

       %Calculate the matrix (n_conduc_param*n_measurements)
       tikhonov_inv= (J'*J +  hp^2*RtR)\J';
       
       %Create an eidors object for iterations
       eidors_obj('set-cache', inv_model, 'mc_2011_one_step_inv', tikhonov_inv);
       eidors_msg('mc_inv_solve: setting cached value', 3);
   end
end
  
end