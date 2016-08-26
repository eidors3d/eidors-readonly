function [ output_args ] = inv_solve_ptensor_lbfgs( imdl, meas_data )
%INV_SOLVE_PTENSOR_LBFGS
%
% L-BFGS inverse solver with a polarization tensor approximation to diag(H)
% for the initial Hessian at each iteration.

% Hard coded number of eigenvectors to store...
mem = 20;
max_its = 50;

g_tol = 1e-6;
r_tol = 1e-6;

% Last m differences between x_i and x_i-1
S = zeros(numel(m_k), mem);

% Last m differences between grad x_i+1 and grad x_i
Y = zeros(numel(m_k), mem);

% Indexing for S/Y
m_ind = repmat(1:mem, 1, ceil((max_its + 2)/2));



% Initialise 2-loop parameter
alpha_i = zeros(1,mem);

% Initialise res vector
resvec = zeros(params.outer.max_its + 1,1);

% Initialise norm of gradient
g = resvec;

d_phi_0 = 0;


% Stop if either gradient or residual reach predefined tol
% Do while (||k==1) to force always one iteration
while ( g(k) > g_tol  && resvec(k)/resvec(1) > r_tol ) || k==1
    fprintf('Beginning iteration %i:\n',k)
    
    % ---------------------------------------------------------------------
    % Implicit Hessian approximation
    disp('  calculating descent direction...')
           
    if k==1
        if isfield(params.outer, 'gamma0')
            gamma_k = params.outer.gamma0;
        else
            gamma_k = 1.;
        end
    else
        gamma_k = ( S(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) )/ ...
            ( Y(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) );
    end

    
    % 2-loop recursions for implicit Hessian approximation
    q = G_new;
    
    % ---------------
    % Loop 1
    for i=k-1:-1:max(k-mem_k,1)
        alpha_i(i) = S(:,m_ind(i)).' * q/( Y(:,m_ind(i)).'*S(:,m_ind(i)) );
        q = q - alpha_i(i) * Y(:,m_ind(i));
    end
    
    % Initial Hessian -----------------------
    H0 = calc_phessian_obj(imdl, img, DU0, DN, delta_d );
    
    % Ensure positivity
    H0(H0<h_min) = h_min;
    
    p_k = 1./H0;
        
       
    
    % ---------------------
    % Loop 2
    for i=max(1,k-mem_k) : k-1
        beta = Y(:,m_ind(i)) .' * p_k/( Y(:,m_ind(i)).' * S(:,m_ind(i)) );
        p_k = p_k + S(:,m_ind(i)) * (alpha_i(i) - beta);
    end
    
        
    % Do linesearch
    
    
    
end



end

