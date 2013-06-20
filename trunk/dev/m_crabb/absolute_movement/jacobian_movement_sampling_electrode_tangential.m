function J = jacobian_movement_sampling_electrode_tangential(fwd_model, img)
% JACOBIAN_MOVEMENT   Computes the Jacobian matrix for conductivity and
% electrode movement variables in 3D EIT.
% Args:     fwd_model - the EIDORS object forward model
%            img - the image background conductivity
%
% fwd_model.conductivity_jacobian - function to calculate conductivity
%                                   Jacobian (defaults to jacobian_adjoint)
%
% Returns:          J - the Jacobian matrix [Jc, Jm]
%
% WARNING: THIS CODE IS EXPERIMENTAL AND GIVES PROBLEMS

if isstr(fwd_model) && strcmp(fwd_model,'UNIT_TEST'); do_unit_test; return ; end

if nargin == 1
   img= fwd_model;
else
   warning('EIDORS:DeprecatedInterface', ...
      ['Calling JACOBIAN_MOVEMENT with two arguments is deprecated and will cause' ...
       ' an error in a future version. First argument ignored.']);
   warning off EIDORS:DeprecatedInterface

end
fwd_model= img.fwd_model;

if isfield(fwd_model,'conductivity_jacobian')
   img.fwd_model.jacobian = fwd_model.conductivity_jacobian;
   Jc= calc_jacobian( img );
else
   img.fwd_model = mdl_normalize(fwd_model, 0); % we normalize on our own
   Jc = jacobian_adjoint(img);
end

Jm = jacobian_electrode_movement_sampling(fwd_model,img);
J=[Jc,Jm];
