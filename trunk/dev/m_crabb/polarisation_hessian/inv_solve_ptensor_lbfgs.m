function [ x_k, resvec, fwd_cts, abs_er ] = inv_solve_ptensor_lbfgs( imdl, img, data, hess_opts, true_im )
%INV_SOLVE_PTENSOR_LBFGS
%
% L-BFGS inverse solver with a polarization tensor approximation to diag(H)
% for the initial Hessian at each iteration.
%
% NB default is no precon

if nargin==3
    hess_opts = [];
end

if isfield(hess_opts,'H0_type')
    H0_type = hess_opts.H0_type;
end
if isfield(hess_opts,'rescale')
    rescale = hess_opts.rescale;
else
    rescale = 0;
end
if isfield(hess_opts,'use_hyper')
    use_hyper = hess_opts.use_hyper;
else
    use_hyper = 0;
end
if isfield(hess_opts, 'update_U0')
    update_U0 = hess_opts.update_U0;
else
    update_U0 = 0;
end
if isfield(hess_opts,'neumann')
    neumann = hess_opts.neumann;
else
    neumann = 'freespace';
end
if isfield(hess_opts,'update_delta')
    update_delta = hess_opts.update_delta;
else
    update_delta = 1;
end

if isfield(hess_opts, 'flexible')
    flexi = hess_opts.flexible;
else
    flexi = true;
end

if isfield(hess_opts,'ptensor_its')
    ptensor_its = hess_opts.ptensor_its;
else
    ptensor_its = inf;
end

if isfield(hess_opts, 'max_its')
    max_its = hess_opts.max_its;
else
    max_its = 50;
end
    
if isfield(hess_opts, 'mem')
    mem = hess_opts.mem;
else
    mem = 20;
end
    

global fwd_cts
fwd_cts = 0;



fmdl = img.fwd_model;

img_0 = img;

% Minimum value on the diagonal (force +ve def)
h_min = 1e-12;

g_tol = 0;%1e-9;
r_tol = 0;%1e-12;
updt_tol = 1e-12;

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

%
if nargin==5
    x_true = true_im.elem_data;
    abs_er(1) = norm(x_true - x_k);
end


% Initialise norm of gradient
g = resvec;

% im_bkgnd = calc_jacobian_bkgnd(imdl);


% Gradient
g_k = calc_grad(x_k, img_0, imdl, data); % this comes out completely wrong!


% J = calc_jacobian(img);
% g_k = J.'*delta_d;

g(1) = norm(g_k);
resvec(1) = objective(x_k, img_0, imdl, data);

k=1;

% Background homogeneous fields
% Switch to update this for new img/do it once?
u0= d_k.volt;
DU0 = calc_grad_potential(img_0, u0); % imdl.img ?

% Stop if either gradient or residual reach predefined tol
% Do while (||k==1) to force always one iteration
while ( g(k) > g_tol  && resvec(k)/resvec(1) > r_tol ) || k==1 % TODO: stop conditions to match EIDORS standards?
    
    
    
    fprintf('Beginning iteration %i:\n',k)
    
    % ---------------------------------------------------------------------
    % Implicit Hessian approximation
    disp('  calculating descent direction...')
         
    
    % 2-loop recursions for implicit Hessia6  n approximation
    q = g_k;
    
    % ---------------
    % Loop 1
    for ii=k-1:-1:max(k-mem,1)
        alpha_i(ii) = S(:,mod(ii-1,mem)+1).' * q/( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1,mem)+1) );
        q = q - alpha_i(ii) * Y(:,mod(ii-1,mem)+1);
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
                            max_C0 = (hp^2*spdiags(RtR) + D0 - h_min)./C0;
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
                % Re-scale
                if k==1
                    gamma_k = 1;
                else
                    gamma_k = ((H0*S(:,mod(k-2,mem)+1)).' * Y(:,mod(k-2,mem)+1) )/ ...
                        ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) );
                end
                
                % Rescale to fit diff in gradients
                if rescale
                    H0 = H0 * gamma_k;
                end
                
            end % Check within update range
            
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
                H0 = calc_hessian_obj(imdl.fwd_model, img, 1:length(x_k), delta_d);
                fwd_cts = fwd_cts + length(x_k) + numel(delta_d);
                
                H0 = spdiags(diag(H0), 0, length(x_k), length(x_k));
                
                if use_hyper
                    RtR = calc_RtR_prior(imdl);
                    hp = calc_hyperparameter(imdl);
                    H0 = H0 + hp^2 * RtR;
                end
                
            end
            
        otherwise
            if k==1
                gamma_k = 1.;
            else
                gamma_k = ( S(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) )/ ...
                    ( Y(:,mod(k-2,mem)+1).' * Y(:,mod(k-2,mem)+1) );
            end
            H0 = max(gamma_k, h_min);
    end % Hess type switch

    
    
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
    if any(isnan(p_k)) || any(isinf(p_k))
        p_k = g_k/gamma_k;
    end
        
    % Do linesearch - this should really satisfy strong Wolfe conditions
    % for BFGS
    
    % Find steplength satisfying strong Wolfe conditions (for +ve definite
    % update)
%     stepl = 1;
    
%     [img, stepl] = line_optimize(img, -p_k, data); % Eidors default
    stepl = wolfe_search(x_k, -p_k, 1, 10, @(x)objective(x, img, imdl, data), @(x)calc_grad(x, img, imdl, data));
    
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
        abs_er(k+1) = norm(x_true - x_k);
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
        fvdsgdf=1; %nonsense for breakpoint
        
    end
    
    g_k = g_k2;
    
    g(k+1) = norm(g_k);
    
    %
    k=k+1;
    
    if norm(p_k*stepl) < updt_tol || k> max_its % k is iter+1
        break
    end
    
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


