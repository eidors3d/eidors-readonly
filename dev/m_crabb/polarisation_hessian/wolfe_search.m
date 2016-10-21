function [ alpha ] = wolfe_search( x, p, alpha_1, alpha_max, F, G )
%WOLFE_SEARCH Summary of this function goes here
%   Wolfe linesearch, including quadratic interpolation



c1 = 1e-4;
c2 = 0.95;
max_it = 50;
alpha_step = 1.25;
delta_min = 1e-6;

% Disable singular matrix warnings -- handled by safeguarding
% interpolations
warning('off','MATLAB:nearlySingularMatrix');

% Initialise
alpha_old = 0; alpha_new = alpha_1; alpha = alpha_1; k=1;

% Descent direction gradient at alpha=0
phi_0 = F(x);
grad_f = G(x);
d_phi_0 = dot(grad_f, p);
phi_old = phi_0;

% Initialise interp pts
Gradpts = zeros(max_it,2);
Gradpts(1,2) = d_phi_0;

% Run optimisation
while k<=max_it
    
    phi_new = F(x + alpha_new * p);    
    
    % New value does not satisfy sufficient decrease condition
    if phi_new> phi_0 + c1*alpha_new* d_phi_0 || (phi_new >= phi_old && k>1)
        alpha = zoom_search(alpha_old, alpha_new, x, p,F,G);
        break
    end % if
    
    % Calc gradient
    grad_f = G( x + alpha_new * p);
    
    % Directional derivative
    d_phi_new =  dot(grad_f,p);
    
    % Satisfies Wolfe conditions - break
    if abs(d_phi_new) <= - c2 * d_phi_0
        alpha = alpha_new;
        break;
    end % if

    % Increasing function
    if d_phi_new >=0
        alpha = zoom_search(alpha_new, alpha_old, x, p, F,G);
        break
    end % if
    

       
    % Save gradient interp pts
    Gradpts(k,1) = alpha_new;
    Gradpts(k,2) = phi_new;
    
    % update alpha
    alpha_old = alpha_new;
        
    % Try interp pt
    if k>=3 
        % Quadratic interp
        alph_interp_grads = quad_interp_vvv(Gradpts(1+k-3,1), Gradpts(2+k-3,1), Gradpts(3+k-3,1),...
                                            Gradpts(1+k-3,2), Gradpts(2+k-3,2), Gradpts(3+k-3,2));
    else
        % Linear interp
        alph_interp_grads = Gradpts(2,2)*Gradpts(2,1)/(Gradpts(1,2) - Gradpts(2,2));

    end
    
      
    % Take the larger step of interpolated and stepped for sufficient progress
    if abs(alph_interp_grads - alpha_old)/alpha_old < delta_min
        alpha_new = max(alpha_step*alpha_new, alph_interp_grads);
    end
    
    % Check for (close to) max step
    if alpha_new > alpha_max || abs(alpha_new - alpha_max)/alpha_max < delta_min
        alpha_new = alpha_max;
    end
    
    

    % update old phi
    phi_old = phi_new;
    
    % update iteration number
    k = k+1;
    
end % while

% Return smallest if none found
if k>max_it
    [~,indx] = min(Gradpts(:,2));
    alpha = Gradpts(indx,1);
end

% re-enable warnings
warning('on','MATLAB:nearlySingularMatrix');


end

% Quadratic interpolation, f(x1), f(x2), f(x3)
function x_min = quad_interp_vvv(x1, x2, x3, v1, v2, v3)
    vals = [x1^2, x1, 1; x2^2, x2, 1; x3^2, x3, 1];
    coefs = vals\[v1; v2; v3];

    % Trial value is minimiser of quadratic
    x_min = -coefs(2)/(2*coefs(1));
end