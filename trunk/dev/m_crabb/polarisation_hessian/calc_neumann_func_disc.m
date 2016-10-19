function [ N ] = calc_neumann_func_disc( x, z )
%CALC_GRAD_N 
% INPUT
% N(x,z) freespace Neumann i.e.
% \nabla \cdot (\nabla N) = \delta
% Single measure loc x, vector of pts z

%Get set of inversion points
z_star = z/(norm(z,2))^2;

% 2D case
if size(z,2)==2 && size(x,2)==2
    
    % Greens function contribution
    RR = x - z;
    RR_star = x-z_star;

    N = - (1/(2*pi))*log(sqrt(sum(RR.^2,2)))...
        - (1/(2*pi))*log(sqrt(sum(RR_star.^2,2)));%...
        %+ (1/(2*pi))*log(norm(z,2));%...
        %+ (1/(4*pi))*(norm(x,2))^2;
        
        
    
    % Correction term on for boundary
    

else
    error('mismatched dimensions for x/z or not in 2D')
end


end

