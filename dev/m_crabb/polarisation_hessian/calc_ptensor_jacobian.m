function [ J ] = calc_ptensor_jacobian( fmdl, img, DU0, type )
%CALC_PTENSOR_JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

%If freespace then do grads for each elements
if nargin==3
    type = 'freespace';
end

% Extract element centres z
zs = fmdl.elem_centre;

% Elem volumes   
vs = pi*fmdl.elem_volume;

ks = img.elem_data;

% Extract tensor components, else assumes tensors are circles
if isfield(fmdl, 'M_tensor')
    a = fmdl.M_tensor.a;
    b = fmdl.M_tensor.b;
    rots = fmdl.M_tensor.rot;
else
    a = ones(size(fmdl.elems,1),1);
    b = ones(size(fmdl.elems,1),1);
    rots = ones(size(fmdl.elems,1),1);
end

% First tensor component (non-differentiated)
A = (a+b)./(a + ks.*b);
B = (a+b)./(b + ks.*a);

P0_11 = A.*cos(rots).^2 - B.*sin(rots).^2; 
P0_12 = (A-B).*cos(rots).*sin(rots);
P0_22 = A.*sin(rots).^2 + B.*cos(rots).^2;

% Second tensor component (first derivative)
A = -(a+b).*b./((a + ks.*b).^2);
B = -(a+b).*a./((b + ks.*a).^2);

P1_11 = A.*cos(rots).^2 - B.*sin(rots).^2;
P1_12 = (A-B).*cos(rots).*sin(rots);
P1_22 = A.*sin(rots).^2 + B.*cos(rots).^2;

% % Second derivative of tensor
% A = 2*(a+b).*b.^2./((a + ks.*b).^3);
% B = 2*(a+b).*a.^2./((b + ks.*a).^3);
% 
% P2_11 = A.*cos(rots).^2 - B.*sin(rots).^2;
% P2_12 = (A-B).*cos(rots).*sin(rots);
% P2_22 = A.*sin(rots).^2 + B.*cos(rots).^2;


% Counter
meas_start = 1;

% Initialise
du2 = zeros(size(zs,1),1);
d2u = du2;

% Loop over source terms
n_drive = length(fmdl.stimulation);

J = cell(n_drive,1);

for ii=1:n_drive

    % Loop over measurements
    M = fmdl.stimulation(ii).meas_pattern;
    [~,meas] = find(M);
    meas = unique(meas);
    
    n_meas = size(M,1);
    
    DU0_M_DN = zeros(size(zs,1),size(M,2));
    DU0_M1_DN = DU0_M_DN;
    DU0_M2_DN = DU0_M_DN;
    
%     J{ii} = zeros(size(meas,1),size(zs));
    
    % NB measurements are not single electrode things...
    for jj=meas.'
        
        % Measurement location - non weighted mean location!!
        meas_loc = fmdl.nodes(fmdl.electrode(jj).nodes,:);
        meas_loc = sum(meas_loc,1)/size(meas_loc,1);
            
        % Tensor shape for elements, or assume circular
%         if isfield(fmdl, 'M_tensor')
            % Only for 2D!   
            % Temp while not being passed
            switch type
                case 'freespace'
                    DN(:,:,ii) = -calc_neumann_grad_func_freespace(meas_loc, zs);
                case 'disc'
                    DN(:,:,ii) = -calc_neumann_grad_func_disc(meas_loc, zs);
                case 'fem_neuman'
            end
            
            % Non-differentiated
            DU0_M_DN(:,jj) = sum(DU0(:,:,ii).*[P0_11.*DN(:,1,ii) + P0_12.*DN(:,2,ii),...
                                         P0_12.*DN(:,1,ii) + P0_22.*DN(:,2,ii)],2);
                                     
            % Second tensor component (first derivative)
            DU0_M1_DN(:,jj) = sum(DU0(:,:,ii).*[P1_11.*DN(:,1,ii) + P1_12.*DN(:,2,ii),...
                                         P1_12.*DN(:,1,ii) + P1_22.*DN(:,2,ii)],2);
            
%             % Second derivative of tensor component            
%             DU0_M2_DN(:,jj) = sum(DU0(:,:,ii).*[P2_11.*DN(:,1,ii) + P2_12.*DN(:,2,ii),...
%                                          P2_12.*DN(:,1,ii) + P2_22.*DN(:,2,ii)],2);

        
    end
    
    % First derivatives
    J{ii} = (repmat(vs, 1, size(M,1)).*(M*(DU0_M_DN.')).'+ repmat(vs.*(ks-1), 1, size(M,1)).*(M*(DU0_M1_DN.')).').';
%     du2 = du2 + sum(drive_contn_D1.^2,2);
    
    
    % Second derivatives
%     drive_contn_D2 = repmat(vs, 1, size(M,1)).*(M*(2*DU0_M1_DN+...
%         repmat(ks-1,1,size(M,2)).*DU0_M2_DN).').'.*repmat(delta_d(meas_start:meas_start+n_meas-1),1,size(zs,1)).';
%     d2u = d2u + sum(drive_contn_D2,2);
    
    meas_start = meas_start + n_meas;

end

J = cell2mat(J);

end

