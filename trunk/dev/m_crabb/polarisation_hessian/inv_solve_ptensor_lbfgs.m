function [ output_args ] = inv_solve_ptensor_lbfgs( imdl, meas_data )
%INV_SOLVE_PTENSOR_LBFGS
%
% L-BFGS inverse solver with a polarization tensor approximation to diag(H)
% for the initial Hessian at each iteration.

% Hard coded number of eigenvectors to store...
mem = 20;
max_its = 50;

% Minimum value on the diagonal (force +ve def)
h_min = 1e-6;

g_tol = 1e-6;
r_tol = 1e-6;

% Background fields for Hess approxn
imdl.fwd_solve.get_all_meas=1;
[d_k] = fwd_solve_higher_order(imdl);
delta_d = d_k.meas - meas_data;
u0= d_k.volt;
DU0 = calc_grad_potential(imdl, u0); % imdl.img ?
H0 = calc_phessian_obj(imdl, img, DU0, DN, delta_d );[d0] = fwd_solve_higher_order(imdl);


u0 = d0.volt;
DN = calc_neuman_grad_func(imdl); 

% Last m differences between x_i and x_i-1
S = zeros(numel(m_k), mem);

% Last m differences between grad x_i+1 and grad x_i
Y = zeros(numel(m_k), mem);

% Indexing for S,Y
m_ind = repmat(1:mem, 1, ceil((max_its + 2)/2));


% Initialise 2-loop parameter
alpha_i = zeros(1,mem);

% Initialise res vector
resvec = zeros(params.outer.max_its + 1,1);

% Initialise norm of gradient
g = resvec;

% -------------------------------------------------------------------------
% This part needs re-writing for EIDORS format/functions

% image data
x_k = imdl.image_data; %?

% Gradient
g_k = grad_obj();% This done with Jacobian elsewhere rather than adjoint?


d_phi_0 = 0;


% Stop if either gradient or residual reach predefined tol
% Do while (||k==1) to force always one iteration
while ( g(k) > g_tol  && resvec(k)/resvec(1) > r_tol ) || k==1 % TODO: stop conditions to match EIDORS standards?
    
    
    
    fprintf('Beginning iteration %i:\n',k)
    
    % ---------------------------------------------------------------------
    % Implicit Hessian approximation
    disp('  calculating descent direction...')
         
    % Is this needed with Hessian approxn? 
    if k==1
        gamma_k = 1.;
    else
        gamma_k = ( S(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) )/ ...
            ( Y(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) );
    end

    
    % 2-loop recursions for implicit Hessian approximation
    q = g_k;
    
    % ---------------
    % Loop 1
    for i=k-1:-1:max(k-mem_k,1)
        alpha_i(i) = S(:,m_ind(i)).' * q/( Y(:,m_ind(i)).'*S(:,m_ind(i)) );
        q = q - alpha_i(i) * Y(:,m_ind(i));
    end
    
    % Initial Hessian -----------------------
    
    
    % Ensure positivity - there are better ways to do this preserving scale
    H0(H0<h_min) = h_min;
    p_k = 1./H0; % * gamma_k ?
    
    % TODO: 
    %   1) needs contribution from the prior
    %   2) if J is already available, could take J^TJ + second derivs from
    %      pHessian only
    
    % ---------------------
    % Loop 2
    for i=max(1,k-mem_k) : k-1
        beta = Y(:,m_ind(i)) .' * p_k/( Y(:,m_ind(i)).' * S(:,m_ind(i)) );
        p_k = p_k + S(:,m_ind(i)) * (alpha_i(i) - beta);
    end
    
        
    % Do linesearch - this should really satisfy strong Wolfe conditions
    % for BFGS
    
    % stepl = lsearch( x_k, -p_k ...)
    
    % Save difference in gradient and img
    S(:, s_ind(k)) = -stepl * p_k;
    
    x_k = x_k - stepl * p_k;    
    g_k2 = grad_obj(); % <------ find EIDORS grad
    Y(:, s_ind(k)) = g_k2 - g_k;
    g_k = g_k2;
    
    
    % Update norm
    [d_k] = fwd_solve_higher_order(imdl);
    delta_d = d_k.meas - meas_data;
    u0= d_k.volt;
    DU0 = calc_grad_potential(imdl, u0); % imdl.img ?
    H0 = calc_phessian_obj(imdl, img, DU0, DN, delta_d );[d0] = fwd_solve_higher_order(imdl);
    resvec(k+1) = norm(delta_d)^2;
    
    %
    k=k+1;
    
end



end

