function Reg= noser_image_prior( inv_model );
% NOSER_IMAGE_PRIOR calculate image prior
% Reg= noser_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% Prior is diag( diag(J'*J)^exponent )
% param is normally 2, this value can be changed by
% setting inv_model.noser_image_prior.exponent=2

% (C) 2005 Andy Adler. Licenced under the GPL Version 2
% $Id: noser_image_prior.m,v 1.1 2006-08-21 17:56:28 aadler Exp $

    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( inv_model.fwd_model, img_bkgnd);

    exponent= 2;
    if isfiels(inv_model,'noser_image_prior');
       exponent= inv_model.noser_image_prior.exponent;
    end

    l_prior= size(J,2);

    % Reg is spdiags(diag(J'*J),0, l_prior, l_prior);
    diag_col sum(J.^exponent,1)';
    Reg = spdiags( diag_col, 0, l_prior, l_prior);

