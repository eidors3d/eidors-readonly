function [ Hii ] = calc_phessian_obj( fmdl, img, DU0, DN, delta_d )

%INPUT
% F(s+ds) = F(s) + DU0 *M *DN


%CALC_PHESSIAN 
%
%


% 




% Extract element centres z
zs = fmdl.elem_centre;

% Elem volumes
vs = fmdl.elem_volume;






% Counter
meas_start = 1;

% Loop over source terms
n_drive = length(fmdl.stimulation);
for ii=1:n_drive

    
    
    
    
    % Loop over measurements
    M = fmdl.stimulation(ii).meas_pattern;
    [~,meas] = find(M);
    meas = unique(meas);
    
    n_meas = size(M,1);
    
    

    
    DU0_M_DN = zeros(size(zs,1),size(M,2));
    
    % NB measurements are not single electrode things...
    for jj=meas.'
        
        % Measurement location - non weighted mean location!!
        meas_loc = fmdl.nodes(fmdl.boundary(jj,:),:);
        meas_loc = sum(meas_loc,1)/size(meas_loc,1);
        
        % calc DN(x, z)
        DU0_M_DN(:,jj) = sum(DU0(:,:,ii).*calc_grad_N(meas_loc, zs),2);
        
        
        
        
    end
    
    % Contributions scaled by delta_d for D2
    coef_D2 = sum((M*DU0_M_DN.').*...
        repmat(delta_d(meas_start:meas_start+n_meas-1),1,size(zs,1)),1).';
    
    % Sum of contributions
    coef_D1 = sum(M*(DU0_M_DN.'),1).';
    
    % Counter
    meas_start = meas_start + n_meas;
    

end

% Volume
coef_D1 = vs.*coef_D1;
coef_D2 = vs.*coef_D2;

% 1st derivatives
ks = img.elem_data;
du =  ( 2./(1+ks) - 2*(ks - 1)./((1+ks).^2)).*coef_D1;

% 2nd derivatives
d2u = (-4./((1+ks).^2) + 4*(ks-1)./((1+ks).^2)).*coef_D2;


% Diagonal of Hessian
Hii = du.^2 + d2u; % NB need to do this


end

