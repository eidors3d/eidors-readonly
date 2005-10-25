function Reg= aa_e_move_image_prior( inv_model );
% AA_E_MOVE_IMAGE_PRIOR calculate image prior
% Reg= aa_e_move_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% Use a diagonal weighting, just that
% [Reg](i,i) = [H'*H](i,i) for H= Jacobian

% $Id: aa_e_move_image_prior.m,v 1.3 2005-10-25 17:53:06 aadler Exp $

fwd_model= inv_model.fwd_model;

homg_img= eidors_obj('image', 'homogeneous image', ...
                     'elem_data', ones( size(fwd_model.elems,1) ,1), ...
                     'fwd_model', fwd_model );

J= calc_jacobian( fwd_model, homg_img );

% We want: [Reg](i,i) = [H'*H](i,i) for H= Jacobian
diag_elems= sum( J.^2 );
idx= 1:length( diag_elems);

Reg= sparse( idx, idx, diag_elems);
