function [RM, J, RtR, W, Jfine] = get_RM( inv_model )
%% obtains the reconstruction matrix from a given inverse model
% Fabian Braun <fbn(ät)csem{döt}ch>, 17/12/2015
if isfield(inv_model, 'solve_use_matrix') && isfield(inv_model.solve_use_matrix, 'RM')
    % GREIT
    RM = inv_model.solve_use_matrix.RM;
else
    % GN
    img_bkgnd= calc_jacobian_bkgnd( inv_model );
    J = calc_jacobian( img_bkgnd);
    if nargout > 4
        if isfield(img_bkgnd.fwd_model, 'coarse2fine')
            img_bkgnd.fwd_model = rmfield(img_bkgnd.fwd_model, 'coarse2fine');
        end
        Jfine = calc_jacobian( img_bkgnd);
    end
    
    RtR = calc_RtR_prior( inv_model );
    W   = calc_meas_icov( inv_model );
    hp  = calc_hyperparameter( inv_model );
    
    RM = (J'*W*J +  hp^2*RtR)\J'*W;
end
end