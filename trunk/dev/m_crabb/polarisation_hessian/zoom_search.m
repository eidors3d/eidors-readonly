function [ alpha ] = zoom_search( alpha_lo, alpha_hi, x, p, F, G  )
%ZOOM_SEARCH Summary of this function goes here
%   Detailed explanation goes here

% Disable singular matrix warnings -- handled by safeguarding
% interpolations
warning('off','MATLAB:nearlySingularMatrix');

c1 = 1e-4;
c2 = 0.95;

phi_0 = F(x);
grad_F = G(x);
d_phi_0 = dot(grad_F, p);

phi_lo = F(x + alpha_lo * p);
phi_hi = F(x + alpha_hi * p);

d_phi_lo = dot(G(x + alpha_lo * p),p);

max_it = 50;
delta_min = 1e-6;

for k =1:max_it
    % Determine range is wide enough for sufficient move /force small step
    % alpha_lo < delta_min corrected outside loop
    if abs(alpha_hi - alpha_lo)/alpha_lo < 0.75*delta_min
        alpha = alpha_lo;
        break
        
    end
    
    % Interpolate to find min using appropriate cubic or quadratic spline
    if k==1
        % Have additional gradient info
        if alpha_lo > 0 && alpha_hi>0
            alpha_int = cubic_interp_vvgg(alpha_lo, alpha_hi, 0, phi_lo, phi_hi, d_phi_lo, d_phi_0);
            
            % Ensure doesn't tend to infinity, interpolation wasn't singular
            % use quadratic if so (further safeguarding follows)
            if alpha_int > max(alpha_lo, alpha_hi) + delta_min*alpha_lo || isnan(alpha_int)
                alpha_int = quad_interp_vvg(alpha_lo, alpha_hi, phi_lo, phi_hi, d_phi_lo);
            end
            
        % No additional gradient info as phi_lo=0
        else
            % Quadratic interpolation minimum
            alpha_int = quad_interp_vvg(alpha_lo, alpha_hi, phi_lo, phi_hi, d_phi_lo);
        end
    else
        % Have phi_old info
        alpha_int = cubic_interp_vvvg(alpha_lo, alpha_hi, alpha_old, phi_lo, phi_hi, phi_old, d_phi_lo);
        
        % Ensure doesn't tend to infinity, interpolation wasn't singular
        % use quadratic if so (further safeguarding follows)
        if alpha_int > max(alpha_lo, alpha_hi) - delta_min*alpha_lo || ...
                alpha_int < min(alpha_lo, alpha_hi) + delta_min*alpha_lo...
                || isnan(alpha_int)
            alpha_int = quad_interp_vvg(alpha_lo, alpha_hi, phi_lo, phi_hi, d_phi_lo);
        end
        
    end
    
    % ---------------------------------------------------------
    % Safeguard to ensure in range, sufficient move, closer to alpha_lo,
    % isnt NaN
    if alpha_int < min(alpha_lo, alpha_hi) + delta_min*alpha_lo || ...
            alpha_int > max(alpha_lo, alpha_hi) - delta_min*alpha_lo ||...
            abs(alpha_int - alpha_lo) > abs(alpha_int - alpha_hi)||...
            isnan(alpha_int)
            
        % Closer to alpha_lo or hi?
        if abs(alpha_lo - alpha_int) < abs(alpha_hi - alpha_int)
            alpha = (3*alpha_lo + alpha_hi)/4;
        else
            alpha = ( alpha_lo + alpha_hi )/2;
        end
        
    % Otherwise acceptable    
    else
        alpha = alpha_int;
        
    end
    
    % -------------------------------------------------------
    % evaluate
    phi_alpha = F(x + alpha *p);
    
    % -------------------------------------------------------
    % Increased, choose step length in (alpha_low, alpha)
    if (phi_alpha > phi_0  + c1 * alpha * d_phi_0 || ...
            phi_alpha >= phi_lo) && k < max_it
        % Save old
        alpha_old = alpha_hi;
        phi_old = phi_hi;
        
        % calculate cost
        alpha_hi = alpha; 
    
   % --------------------------------------------------------
   % Not increasing/last it, check gradient condition
    else
        grad_f = G( x + alpha * p);
        d_phi_alpha = dot(grad_f , p);        
        % Satisfies Wolfe conditions, return value
        if (abs(d_phi_alpha) <= -c2 * d_phi_0)
            break;
        end
        

        % Choose new range (alpha_lo, alpha_hi)
        % Determine which side of alpha min lies
        if d_phi_alpha * (alpha_hi - alpha_lo ) >=0
            % Save old -- alpha_lo being dropped
            alpha_old = alpha_hi;
            phi_old = phi_hi;
            
            % update alpha_hi
            alpha_hi = alpha_lo;
            
        else
            % Save old -- alpha_hi being dropped
            alpha_old = alpha_lo;
            phi_old = phi_lo;
            
        end
        %
        phi_lo = phi_alpha; 
        d_phi_lo = d_phi_alpha;
        alpha_lo = alpha;
    end % if

end

% re-enable warnings
warning('on','MATLAB:nearlySingularMatrix');

end

% Interpolation minima
% Quadratic, f(x1), f(x2), f'(x1)
function x_min = quad_interp_vvg(x1, x2, v1, v2, g1)
vals = [x1^2, x1, 1; x2^2, x2, 1; 2*x1, 1, 0];
    coefs = vals\[v1; v2; g1];
    
    % Trial value is minimiser of quadratic
    x_min = -coefs(2)/(2*coefs(1));
end
 

% Cubic, f(x1), f(x2), f'(x1), f'(x3)
function x_min = cubic_interp_vvgg(x1, x2, x3, v1, v2, g1, g3)
vals = [x1^3, x1^2, x1, 1;...
        x2^3, x2^2, x2, 1;...
        3*x1^2, 2*x1, 1, 0;...
        3*x3^2, 2*x3, 1, 0];
    coefs = vals\[v1; v2; g1; g3];
    
    % Trial value is minimiser of quadratic
    x_min = (-coefs(2) + sqrt(coefs(2)^2 - 3*coefs(1)*coefs(3)))/(3*coefs(1));
end

% Cubic, f(x1), f(x2), f(x3), f'(x1)
function x_min = cubic_interp_vvvg(x1, x2, x3, v1, v2, v3, g1)
vals = [x1^3, x1^2, x1, 1;...
        x2^3, x2^2, x2, 1;...
        x3^3, x3^2, x3, 1;...
        3*x1^2, 2*x1, 1, 0];
    coefs = vals\[v1; v2; v3; g1];
    
    % Trial value is minimiser of quadratic
    x_min = (-coefs(2) + sqrt(coefs(2)^2 - 3*coefs(1)*coefs(3)))/(3*coefs(1));
end
