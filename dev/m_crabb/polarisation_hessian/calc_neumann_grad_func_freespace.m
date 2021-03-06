function [ DN ] = calc_neumann_grad_func_freespace( x, z )
%CALC_GRAD_N 
% INPUT
% Gradient of N(x,z) wrt z
% Single measure loc x, vector of pts z


% 2D case
if size(z,2)==2 && size(x,2)==2
    
    % Greens function contribution
    RR = repmat(x, size(z,1),1) - z;
   % DN = RR./(2*pi*repmat(sqrt(sum(RR.^2,2)),1,2));
    DN = RR./(2*pi*repmat(sum(RR.^2,2),1,2));
    
    % Correction term on for boundary
    

else
    error('mismatched dimensions for x/z or not in 2D')
end


end

