function [ Hii ] = calc_phessian_obj( fmdl, img, DU0, DN, delta_d )

%INPUT
% F(s+ds) = F(s) + DU0 *M *DN


%CALC_PHESSIAN 
%
%


% 

% Temp while DN function not complete
DN = zeros(size(DU0));

% Extract element centres z
zs = fmdl.elem_centre;

% Elem volumes
vs = fmdl.elem_volume;

ks = img.elem_data;

% Extract tensor components, else assumes tensors are circles
if isfield(fmdl, 'M_tensor')
    a = fmdl.M_tensor.a;
    b = fmdl.M_tensor.b;
    rots = fmdl.M_tensor.rot;
    
    % First tensor component (non-differentiated)
    A = (a+b)./(a + ks.*b);
    B = (a+b)./(b + ks.*a);
    
    P0_11 = A.*cos(rots).^2 - B.*sin(rots).^2;
    P0_12 = (A-B).*cos(rots).*sin(rots);
    P0_22 = A.*sin(rots).^2 - B.*cos(rots).^2;
    
     % Second tensor component (first derivative)
     A = -(a+b).*b./((a + ks.*b).^2);
     B = -(a+b).*a./((a + ks.*b).^2);
     
     P1_11 = A.*cos(rots).^2 - B.*sin(rots).^2;
     P1_12 = (A-B).*cos(rots).*sin(rots);
     P1_22 = A.*sin(rots).^2 - B.*cos(rots).^2;
     
     % Second derivative of tensor
     A = (a+b).*b.^2./((a + ks.*b).^3);
     B = (a+b).*a.^2./((a + ks.*b).^2);
     
     P2_11 = A.*cos(rots).^2 - B.*sin(rots).^2;
     P2_12 = (A-B).*cos(rots).*sin(rots);
     P2_22 = A.*sin(rots).^2 - B.*cos(rots).^2;
    
end


% Counter
meas_start = 1;

% Initialise
coef_D1 = zeros(size(zs,1),1);
coef_D2 = coef_D1;

% Loop over source terms
n_drive = length(fmdl.stimulation);
for ii=1:n_drive

    % Loop over measurements
    M = fmdl.stimulation(ii).meas_pattern;
    [~,meas] = find(M);
    meas = unique(meas);
    
    n_meas = size(M,1);
    
    DU0_M_DN = zeros(size(zs,1),size(M,2));
    DU0_M1_DN = DU0_M_DN;
    DU0_M2_DN = DU0_M_DN;
    
    
    
    % NB measurements are not single electrode things...
    for jj=meas.'
        
        % Measurement location - non weighted mean location!!
        meas_loc = fmdl.nodes(fmdl.boundary(jj,:),:);
        meas_loc = sum(meas_loc,1)/size(meas_loc,1);
            
        % Tensor shape for elements, or assume circular
        if isfield(fmdl, 'M_tensor')
            % Only for 2D!
   
            % Temp while not being passed
            DN(:,:,ii) = calc_neumann_grad_func_freespace(meas_loc, zs);
            
            % Non-differentiated
            DU0_M_DN(:,jj) = sum(DU0(:,:,ii).*[P0_11.*DN(:,1,ii) + P0_12.*DN(:,2,ii),...
                                         P0_12.*DN(:,1,ii) + P0_22.*DN(:,2,ii)],2);
                                     
            % Second tensor component (first derivative)
            DU0_M1_DN(:,jj) = sum(DU0(:,:,ii).*[P1_11.*DN(:,1,ii) + P1_12.*DN(:,2,ii),...
                                         P1_12.*DN(:,1,ii) + P1_22.*DN(:,2,ii)],2);
            
            % Second derivative of tensor component            
            DU0_M2_DN(:,jj) = sum(DU0(:,:,ii).*[P2_11.*DN(:,1,ii) + P2_12.*DN(:,2,ii),...
                                         P2_12.*DN(:,1,ii) + P2_22.*DN(:,2,ii)],2);
                                     
        else
            
            % calc DN(x, z)
            DU0_M_DN(:,jj) = sum(DU0(:,:,ii).*calc_neumann_grad_func_freespace(meas_loc, zs),2);
        
        end
        
        
    end
    
    
    % TODO: add contributions from each ii, currently only saving last?
    
    if isfield(fmdl, 'M_tensor')
        % First derivatives
        coef_D1 = coef_D1 + sum(M*(DU0_M_DN.'),1).' + (ks-1).*(sum(M*(DU0_M1_DN.'),1).');
               
        % Second derivatives
        coef_D2 = coef_D2 + sum(-2*M*(DU0_M1_DN.').*repmat(delta_d(meas_start:meas_start+n_meas-1),1,size(zs,1)),1).'...
            + (ks - 1).*sum(M*(DU0_M2_DN.').*repmat(delta_d(meas_start:meas_start+n_meas-1),1,size(zs,1)),1).';
        
        
    else
        
        % Contributions scaled by delta_d for D2
        coef_D2 = coef_D2 + sum((M*DU0_M_DN.').*...
            repmat(delta_d(meas_start:meas_start+n_meas-1),1,size(zs,1)),1).';
        
        % Sum of contributions
        coef_D1 = coef_D1 + sum(M*(DU0_M_DN.'),1).';
        
       
    end
    
     % Counter
        meas_start = meas_start + n_meas;

end

% Volume
coef_D1 = vs.*coef_D1;
coef_D2 = vs.*coef_D2;

if isfield(fmdl, 'M_tensor')
    % Already multiplied with (k-1) factors
    du = coef_D1;
    d2u = coef_D2;
    
else 
    % 1st derivatives
    
    du =  ( 2./(1+ks) - 2*(ks - 1)./((1+ks).^2)).*coef_D1;
    
    % 2nd derivatives
    d2u = (-4./((1+ks).^2) + 4*(ks-1)./((1+ks).^2)).*coef_D2;

end

% Diagonal of Hessian
Hii = du.^2 + d2u; % NB need to do this


end

