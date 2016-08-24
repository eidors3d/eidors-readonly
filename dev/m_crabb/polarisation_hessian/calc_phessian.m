function [ Hii ] = calc_phessian( fmdl, img )
%CALC_PHESSIAN 
%
%


% 




% Extract element centres z
zs = fmdl.elem_centre;

% Elem volumes
vs = fmdl.elem_volume;




% Coeffs
coef_H = zeros(length(zs),1);

% Loop over source terms
n_drive = length(fmdl.stimulation);
for ii=1:n_drive
    % Gradient of Homogeneous solution
    DU0 = zeros(size(zs)); % Not a solution!!
    
    
    
    
    % Loop over measurements
    n_meas = size(fmdl.stimulation(ii).meas_pattern,1);
    
    % NB measurements are not single electrode things...
    for jj=1:n_meas
        
        % calc DN(x, z)
        DN = calc_grad_N(x(jj,:), zs); % Need some x!!
        
        % Contribution
        coef_H = coef_H + sum(DN.*DU0,2);
        
    end

end

% Volume
coef_H = vs.*coef_H;

% 1st derivatives
ks = img.elem_data;
du =  ( 2./(1+ks) - 2*(ks - 1)./((1+ks).^2)).*coef_H;

% 2nd derivatives
d2u = (-4./((1+ks).^2) + 4*(ks-1)./((1+ks).^2)).*coef_H;


% Diagonal of Hessian
Hii = du.^2 + d2u.*(d_sim - d_obs); % NB need to do this


end

