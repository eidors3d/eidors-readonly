function [ output_args ] = inv_solve_ptensor_lbfgs( imdl, img, data )
%INV_SOLVE_PTENSOR_LBFGS
%
% L-BFGS inverse solver with a polarization tensor approximation to diag(H)
% for the initial Hessian at each iteration.

% Hard coded number of eigenvectors to store...
mem = 20;
max_its = 50;


% Minimum value on the diagonal (force +ve def)
h_min = 1e-12;

g_tol = 1e-12;
r_tol = 1e-12;

meas_data = data.meas;

% Background fields for Hess approxn
imdl.fwd_solve.get_all_meas=1;
[d_k] = fwd_solve_higher_order(img);
delta_d = d_k.meas - meas_data;
[d0] = fwd_solve_higher_order(img);


u0 = d0.volt;
% DN = calc_neuman_grad_func(img); 
% image data
x_k = img.elem_data; %?

% Last m differences between x_i and x_i-1
S = zeros(length(x_k), mem);

% Last m differences between grad x_i+1 and grad x_i
Y = zeros(length(x_k), mem);




% Initialise 2-loop parameter
alpha_i = zeros(1,mem);

% Initialise res vector
resvec = zeros(max_its + 1,1);

% Initialise norm of gradient
g = resvec;

% -------------------------------------------------------------------------
% This part needs re-writing for EIDORS format/functions



% Gradient <- change grad to adjoint method?6
J = calc_jacobian(img);
g_k = J.'*delta_d;

g(1) = norm(g_k);
resvec(1) = norm(delta_d);

k=1;

% Stop if either gradient or residual reach predefined tol
% Do while (||k==1) to force always one iteration
while ( g(k) > g_tol  && resvec(k)/resvec(1) > r_tol ) || k==1 % TODO: stop conditions to match EIDORS standards?
    
    
    
    fprintf('Beginning iteration %i:\n',k)
    
    % ---------------------------------------------------------------------
    % Implicit Hessian approximation
    disp('  calculating descent direction...')
         
%     % Is this needed with Hessian approxn? 
%     if k==1
%         gamma_k = 1.;
%     else
%         gamma_k = ( S(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) )/ ...
%             ( Y(:,m_ind(k-1)).' * Y(:,m_ind(k-1)) );
%     end

    
    % 2-loop recursions for implicit Hessia6  n approximation
    q = g_k;
    
    % ---------------
    % Loop 1
    for ii=k-1:-1:max(k-mem,1)
        alpha_i(ii) = S(:,mod(ii-1,mem)+1).' * q/( Y(:,mod(ii-1,mem)+1).'*S(:,mod(ii-1,mem)+1) );
        q = q - alpha_i(ii) * Y(:,mod(ii-1,mem)+1);
    end
    
    % Initial Hessian -----------------------
    % Update phess
    u0= d_k.volt;
    DU0 = calc_grad_potential(img, u0); % imdl.img ?
    H0 = calc_phessian_obj(imdl, img, DU0, [], delta_d, 1 );
    
    % Ensure positivity - there are better ways to do this preserving scale
    H0(H0<h_min) = h_min; % Should this be rescaled also?
    p_k = q./H0; 
    
    % TODO: 
    %   1) needs contribution from the prior
    %   2) if J is already available, could take J^TJ + second derivs from
    %      pHessian only
    
    % ---------------------
    % Loop 2
    for ii=max(1,k-mem) : k-1
        beta = Y(:,mod(ii-1,mem)+1) .' * p_k/( Y(:,mod(ii-1,mem)+1).' * S(:,mod(ii,mem)) );
        p_k = p_k + S(:,mod(ii-1, mem)+1) * (alpha_i(ii) - beta);
    end
    
        
    % Do linesearch - this should really satisfy strong Wolfe conditions
    % for BFGS
    
    % stepl = lsearch( x_k, -p_k ...)
    stepl = 1;
    
    % Save difference in gradient and img
    S(:, mod(k-1,mem)+1) = -stepl * p_k;
    
    % Update
    x_k = x_k - stepl * p_k; 
    img.elem_data = x_k;
    
    % Forward solve
    d_k = fwd_solve_higher_order(img);
    delta_d = d_k.meas - meas_data;
    resvec(k+1) = norm(delta_d);
    
    % New grad
    J = calc_jacobian(img);
    g_k2 = J.'*delta_d;
    Y(:, mod(k-1,mem)+1) = g_k2 - g_k;
    g_k = g_k2;
    g(k+1) = norm(g_k);
    
    
    %
    k=k+1;
    
end



end

