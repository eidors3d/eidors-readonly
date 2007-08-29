function Reg= noser_image_prior( inv_model );
% NOSER_IMAGE_PRIOR calculate image prior
% Reg= noser_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% Prior is diag( diag(J'*J)^exponent )
% param is normally .5, this value can be changed by
% setting inv_model.noser_image_prior.exponent= new_value

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: noser_image_prior.m,v 1.3 2007-08-29 09:11:47 aadler Exp $

    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( inv_model.fwd_model, img_bkgnd);

    exponent= 0.5;
    if isfield(inv_model,'noser_image_prior');
       exponent= inv_model.noser_image_prior.exponent;
    end

    l_prior= size(J,2);

    % Reg is spdiags(diag(J'*J),0, l_prior, l_prior);
    diag_col= sum(J.^2,1)';
    Reg = spdiags( diag_col.^exponent, 0, l_prior, l_prior);

