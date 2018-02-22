function [ img, resvec, fwd_cts, abs_er, H_BFGS ] = inv_solve_ptensor_lbfgs( imdl, img, data, opt, true_im )
%INV_SOLVE_PTENSOR_LBFGS
%
% L-BFGS inverse solver with a polarization tensor approximation to diag(H)
% for the initial Hessian at each iteration.
%
% NB default is no precon

if nargin==3
    opt = [];
end

if isfield(opt,'H0_type')
    H0_type = opt.H0_type;
else
    H0_type = 'ptensor';
end
if isfield(opt,'rescale')
    rescale = opt.rescale;
else
    rescale = 0;
end
if isfield(opt,'use_hyper')
    use_hyper = opt.use_hyper;
else
    use_hyper = true;
end
if isfield(opt, 'update_U0')
    update_U0 = opt.update_U0;
else
    update_U0 = 0;
end
if isfield(opt,'neumann')
    neumann = opt.neumann;
else
    neumann = 'freespace';
end
if isfield(opt,'update_delta')
    update_delta = opt.update_delta;
else
    update_delta = 1;
end

if isfield(opt, 'flexible')
    flexi = opt.flexible;
else
    flexi = true;
end

if isfield(opt,'ptensor_its')
    ptensor_its = opt.ptensor_its;
else
    ptensor_its = inf;
end

if isfield(opt, 'max_its')
    max_its = opt.max_its;
else
    max_its = 50;
end
    
if isfield(opt, 'mem')
    mem = opt.mem;
else
    mem = 20;
end
    

if isfield(opt,'inv_crime')
    inv_crime = opt.inv_crime;
else
    inv_crime =0;
end

global fwd_cts
fwd_cts = 0;



fmdl = img.fwd_model;

img_0 = img;

% Minimum value on the diagonal (force +ve def)
h_min = 1e-8;
try h_min = opt.h_min; end

g_tol = 0;%1e-9;
d_tol = 1e-4; % Same as EIDORS default check
try g_tol = opt.g_tol; end
try d_tol = opt.d_tol; end

updt_tol = 0;%1e-14;

meas_data = data.meas;

% Background fields for Hess approxn
img.fwd_solve.get_all_meas=1;
[d_k] = fwd_solve(img);
delta_d = d_k.meas - meas_data;


% DN = calc_neuman_grad_func(img); 
% image data
x_k = img.elem_data;

% Last m differences between x_i and x_i-1im
S = zeros(length(x_k), mem);

% Last m differences between grad x_i+1 and grad x_i
Y = zeros(length(x_k), mem);




% Initialise 2-loop parameter
alpha_i = zeros(1,mem);

% Initialise res vector
resvec = zeros(max_its + 1,1);
abs_er = resvec;

if(inv_crime==1) && nargin==5
    if nargin==5
        x_true = true_im.elem_data;
        abs_er(1) = norm(x_true - x_k);
    end
else
    c2f = mk_coarse_fine_mapping(true_im.fwd_model,img.fwd_model);
    %
    if nargin==5
        x_true = true_im.elem_data;
        abs_er(1) = norm(x_true - c2f*x_k);
    end
end



% Initialise norm of gradient
g = resvec;

% im_bkgnd = calc_jacobian_bkgnd(imdl);


% Gradient
g_k = calc_grad(x_k, img_0, imdl, data);

% J = calc_jacobian(img);
% g_k = J.'*delta_d;

g(1) = norm(g_k);
resvec(1) = objective(x_k, img_0, imdl, data);

k=1;
rel_red = 0;

% Background homogeneous fields
% Switch to update this for new img/do it once?
u0= d_k.volt;
DU0 = calc_grad_potential(img_0, u0); % imdl.img ?

% Stop if either gradient or residual reach predefined tol
% Do while (||k==1) to force always one iteration
while ( g(k) > g_tol  && rel_red > d_tol ) || k==1 || k==2 % TODO: stop conditions to match EIDORS standards?
    
    
    
    fprintf('Beginning iteration %i:\n',k)
    
    % ---------------------------------------------------------------------
    % Implicit Hessian approximation
    disp('  calculating descent direction...')
         
    
    % 2-loop recursions for implicit Hessian approximation
    q = g_k;
    
    % ---------------
    % Loop 1
    for ii=k-1:-1:max(k-mem,1)
        alpha_i(ii) = S(:,mod(ii-1,mem)+1).' * q/( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1,mem)+1) );
        q = q - alpha_i(ii) * Y(:,mod(ii-1,mem)+1);
    end
    
    if k==1
        gamma_k = 1;
    else
        gamma_k = ((H0*S(:,mod(k-2,mem)+1)).' * Y(:,mod(k-2,mem)+1) )/ ...
            ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) );
    end
    
    % Initial Hessian -----------------------
    switch H0_type
        case 'ptensor'
            if  k<= ptensor_its  && (flexi || k==1)
                
                if update_U0
                    d_k = fwd_solve(img); % doesn't contribute to cts as cached value here
                    u0= d_k.volt;
                    DU0 = calc_grad_potential(img, u0);
                end
                
                [H0, D0, C0] = calc_phessian_obj(fmdl, img_0, DU0, delta_d, neumann );
                
                H0 = spdiags(H0, 0, length(H0), length(H0));% + hp^2 * RtR;
                % +ve def checks should include RtR
                
                % include regularisation
                if use_hyper
                    RtR = calc_RtR_prior(imdl);
                    hp = calc_hyperparameter(imdl);
                    H0 = H0 + hp^2*RtR;
                    
                    % Ensure suff +ve def
                    if any(spdiags(H0,0) < h_min)
                        % Add only amount of 2nd derivatives retaining +ve def if possible
                        if all(D0 + hp^2*spdiags(RtR,0)>h_min)
                            max_C0 = (hp^2*spdiags(RtR,0) + D0 - h_min)./C0;
                            max_C0(isinf(max_C0) | isnan(max_C0)) = inf;
                            max_C0(max_C0<0) = 0;
                            max_C0 = min([max_C0;1]);
                            H0 = spdiags(D0+ max_C0 * C0, 0, length(H0), length(H0)) + hp^2*RtR;
                            
                        else
                            % Gauss Newton part too small
                            indx = D0 + hp^2 * spdiags(RtR,0) < h_min;
                            H0 = spdiags(D0, 0, length(H0), length(H0)) + hp^2*RtR;
                            H0(indx,indx) = h_min;
                            
                        end
                    end
                    
                else
                    
                    % Ensure suff +ve def
                    if any(H0 < h_min)
                        % Add only amount of 2nd derivatives retaining +ve def if possible
                        if all(D0>h_min)
                            max_C0 = (D0 - h_min)./C0;
                            max_C0(isinf(max_C0) | max_C0 < 0 | isnan(max_C0)) = inf;
                            max_C0 = min(max_C0);
                            H0 = D0 + max_C0 * C0;
                            
                        else
                            % Gauss Newton part too small
                            H0 = D0;
                            H0(H0< h_min) = h_min;
                            spdiags(H0, 0, length(H0), length(H0));
                        end
                    end
                end
                
                
                % Rescale to fit diff in gradients
                if rescale
                    H0 = H0 * gamma_k;
                end
                
            end % Check within update range
            
            
        case 'ptensor_full'
            if  k<= ptensor_its  && (flexi || k==1)
                
                if update_U0
                    d_k = fwd_solve(img); % doesn't contribute to cts as cached value here
                    u0= d_k.volt;
                    DU0 = calc_grad_potential(img, u0);
                end
                
                [~, ~, C0, J0] = calc_phessian_obj(fmdl, img_0, DU0, delta_d, neumann );
                
                H0 = spdiags(C0, 0, length(C0), length(C0)) + J0.'*J0;                
                % +ve def checks should include RtR
                
                % include regularisation
                if use_hyper
                    RtR = calc_RtR_prior(imdl);
                    hp = calc_hyperparameter(imdl);
                    H0 = H0 + hp^2*RtR;
                    
                    % Ensure suff +ve def
                    if any(spdiags(H0,0) < h_min)
                        % Add only amount of 2nd derivatives retaining +ve def if possible
                        if all(D0 + hp^2*spdiags(RtR,0)>h_min)
                            max_C0 = (hp^2*spdiags(RtR,0) + D0 - h_min)./C0;
                            max_C0(isinf(max_C0) | isnan(max_C0)) = inf;
                            max_C0(max_C0<0) = 0;
                            max_C0 = min([max_C0;1]);
                            H0 = spdiags(D0+ max_C0 * C0, 0, length(H0), length(H0)) + hp^2*RtR;
                            
                        else
                            % Gauss Newton part too small
                            indx = D0 + hp^2 * spdiags(RtR,0) < h_min;
                            H0 = spdiags(D0, 0, length(H0), length(H0)) + hp^2*RtR;
                            H0(indx,indx) = h_min;
                            
                        end
                    end
                    
                else
                    
                    % Ensure suff +ve def
                    if any(H0 < h_min)
                        % Add only amount of 2nd derivatives retaining +ve def if possible
                        if all(D0>h_min)
                            max_C0 = (D0 - h_min)./C0;
                            max_C0(isinf(max_C0) | max_C0 < 0 | isnan(max_C0)) = inf;
                            max_C0 = min(max_C0);
                            H0 = D0 + max_C0 * C0;
                            
                        else
                            % Gauss Newton part too small
                            H0 = D0;
                            H0(H0< h_min) = h_min;
                            
                        end
                    end
                end
            end
        case 'DGN0'
            % diag Gauss Newton
            if flexi || k==1
                J = calc_jacobian(img);
                fwd_cts = fwd_cts + numel(delta_d);
                
                H0_DGN = sum(J.*conj(J),1);
                H0_DGN = spdiags(H0_DGN.', 0, length(H0_DGN), length(H0_DGN));
                H0 = H0_DGN;
                
                if use_hyper
                    RtR = calc_RtR_prior(imdl);
                    hp = calc_hyperparameter(imdl);
                    H0 = H0 + hp^2 * RtR;
                end
            end
            
        case 'true'
            % diagonal of true Hessian
            if flexi || k==1
                H0 = calc_hessian_diag_obj(imdl.fwd_model, img, 1:length(x_k), delta_d);
                fwd_cts = fwd_cts + length(x_k) + numel(delta_d);
                
                H0 = spdiags(diag(H0), 0, length(x_k), length(x_k));
                
                if use_hyper
                    RtR = calc_RtR_prior(imdl);
                    hp = calc_hyperparameter(imdl);
                    H0 = H0 + hp^2 * RtR;
                end
                
            end
            
        case 'identity'
%             if k==1
%                 gamma_k = 1.;
%             else
%                 gamma_k = ( S(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) )/ ...
%                     ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) );
%             end
            H0 = max(gamma_k, h_min);
            
        case 'regu'
            if k==1
                gamma_k = 1.;
            else
                gamma_k = ( S(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) )/ ...
                    ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) );
            end
            
            % RtR + some identity (RtR not pos def for Laplace)
            H0 = calc_hyperparameter(imdl)^2*calc_RtR_prior(imdl);
            H0 = max(gamma_k, h_min) * eye(size(H0));
            
    end % Hess type switch

    if any(size(H0)==1)
        H0 = spdiags(H0, 0, length(H0), length(H0));
    end
    
    p_k = H0\q; 
    
    % TODO: 
    %   1) needs contribution from the prio
    %   2) if J is already available, could take J^TJ + second derivs from
    %      pHessian only
    
    % ---------------------
    % Loop 2
    for ii=max(1,k-mem) : k-1
        beta = Y(:,mod(ii-1,mem)+1) .' * p_k/( Y(:,mod(ii-1,mem)+1).' * S(:,mod(ii-1,mem)+1));
        p_k = p_k + S(:,mod(ii-1, mem)+1) * (alpha_i(ii) - beta);
    end
    
    % Safeguard
    if any(isnan(p_k)) || any(isinf(p_k)) || dot(p_k, g_k) <0
        warning('lbfgs_solve: BFGS not defined')
        p_k = g_k/gamma_k;
    end
        
    % Do linesearch - this should really satisfy strong Wolfe conditions
    % for BFGS
    
    % Find steplength satisfying strong Wolfe conditions (for +ve definite
    % update)
%     stepl = 1;
    
%     [img, stepl] = line_optimize(img, -p_k, data); % Eidors default
    stepl = wolfe_search(x_k, -p_k, 1, 10, @(x)objective(x, img, imdl, data), @(x)calc_grad(x, img, imdl, data), opt);
    
    % Save difference in gradient and img
    S(:, mod(k-1,mem)+1) = -stepl * p_k;
    
    % Update
    x_k = x_k - stepl * p_k; 
    img.elem_data = x_k;
    
    % Forward solve
%     d_k = fwd_solve(img);
%     delta_d = d_k.meas - meas_data;
    if update_delta
        [resvec(k+1), delta_d] = objective(x_k, img, imdl, data);
    else
        resvec(k+1) = objective(x_k, img, imdl, data);
    end
    
    if nargin==5
        if(inv_crime==1)        
            abs_er(k+1) = norm(x_true - x_k);
        else
            abs_er(k+1) = norm(x_true - c2f*x_k);            
        end
    end
    
    if resvec(k+1) > resvec(k)
        x_k = x_k + stepl * p_k;
        break
    end
    
    % Save
    img.elem_data = x_k;
    
    % New grad
    g_k2 = calc_grad(x_k, img, imdl, data);
    
    
    % check satisfied (some) Wolfe condition
    d_phi_0 = dot(g_k, -p_k);
    d_phi_new = dot(g_k2, -p_k);
    if abs(d_phi_new) < - d_phi_0
    
        Y(:, mod(k-1,mem)+1) = g_k2 - g_k;
        
    else
        warning('lbfgs_solve: stepl did not satisfy Wolfe condition')
        break
        
    end
    
    g_k = g_k2;
    
    g(k+1) = norm(g_k);
    
    %
    k=k+1;
    
    rel_red = (resvec(k-1)-resvec(k))/resvec(1);
    if norm(p_k*stepl) < updt_tol || k> max_its % k is iter+1
        break
    end
    
end


img.elem_data = x_k;
img.current_params = 'conductivity';
resvec = resvec(1:k);


if nargout>=5
    % Use whatever version Hessian used
    H_BFGS = full(H0);
    
%     % Compare with identity
%     H_IBFGS = ( S(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) )/ ...
%                     ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) ) * eye(size(H_BFGS));
%                 
%     H_RBFGS = calc_RtR_prior(imdl) * calc_hyperparameter(imdl)^2;

    for ii=max(1,k-mem) : k-1
        H_BFGS = H_BFGS +  Y(:,mod(ii-1,mem)+1)* Y(:,mod(ii-1,mem)+1).' / ...
            ( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1, mem)+1)) - ...
            H_BFGS * S(:,mod(ii-1, mem)+1) * S(:,mod(ii-1, mem)+1).' * H_BFGS.' / ...
            (S(:,mod(ii-1, mem)+1).'* H_BFGS * S(:,mod(ii-1, mem)+1));
        
        
%         H_IBFGS = H_IBFGS +  Y(:,mod(ii-1,mem)+1)* Y(:,mod(ii-1,mem)+1).' / ...
%             ( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1, mem)+1)) - ...
%             H_IBFGS * S(:,mod(ii-1, mem)+1) * S(:,mod(ii-1, mem)+1).' * H_IBFGS.' / ...
%             (S(:,mod(ii-1, mem)+1).'* H_IBFGS * S(:,mod(ii-1, mem)+1));
%         
%         
%         H_RBFGS = H_RBFGS +  Y(:,mod(ii-1,mem)+1)* Y(:,mod(ii-1,mem)+1).' / ...
%             ( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1, mem)+1)) - ...
%             H_RBFGS * S(:,mod(ii-1, mem)+1) * S(:,mod(ii-1, mem)+1).' * H_RBFGS.' / ...
%             (S(:,mod(ii-1, mem)+1).'* H_RBFGS * S(:,mod(ii-1, mem)+1));
        
        
        
    end
    
%     H_vers = {H_BFGS, H_IBFGS, H_RBFGS};
    
end


end



% Cost
function [C, delta_d] = objective(x_k, images, invmdl, d)
    global fwd_cts
    fwd_cts = fwd_cts + 1;

    im = images;
    im.elem_data = x_k;
    d_im = fwd_solve(im);
    RtR = calc_RtR_prior(invmdl);
    hp = calc_hyperparameter(invmdl);
    imdif = im.elem_data - images.elem_data;
    delta_d = d_im.meas - d.meas;
    C = 0.5*norm(delta_d)^2 + hp^2 * imdif.'*RtR*imdif; 
    

end

% Grad
function G = calc_grad(x_k, image, invmdl, d)
    global fwd_cts
    fwd_cts = fwd_cts + 1;

    im = image;
    im.elem_data = x_k;
    d_im = fwd_solve(im);
    J = calc_jacobian(im);
    RtR = calc_RtR_prior(invmdl);
    hp = calc_hyperparameter(invmdl);
    imdif = x_k - image.elem_data;
    delta_d = d_im.meas - d.meas;
    G = J.'*(delta_d) + hp^2 * RtR*(imdif); % need to generalise prior
end


