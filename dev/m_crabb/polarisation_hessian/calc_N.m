function [ N ] = calc_N( x, z )
%CALC_GRAD_N 
% INPUT
% N(x,z) freespace Neumann
% Single measure loc x, vector of pts z


% 2D case
if size(z,2)==2 && size(x,2)==2
    
    % Greens function contribution
    RR = x - z;
    N = 1/(2*pi*log(sqrt(sum(RR.^2,2))));
    
    % Correction term on for boundary
    

else
    error('mismatched dimensions for x/z or not in 2D')
end


end

