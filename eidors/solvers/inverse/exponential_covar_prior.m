function Reg= exponential_covar_prior( inv_model );
% EXPONENTIAL_COVAR_PRIOR image prior with exponential
%      inter-element covariance, multiplied by a NOSER
%      type weighting, if desired.
% Reg= exponential_covar_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters: exponential rate
%   gamma= inv_model.exponential_covar_prior.gamma
% Parameters: NOSER exponent: diag(J'*J)^exponent )
%   gamma= inv_model.exponential_covar_prior.noser_exp
%       DEFAULT is 0 (ie. no NOSER exponent)
%
% FIXME: this code does an inv or the reg matrix. This
%        is very inefficient.

% (C) 2007 Andy Adler. License: GPL version 2 or version 3
% $Id$

fwd_model= inv_model.fwd_model;

try 
    gamma= inv_model.exponential_covar_prior.gamma;
catch
    xy_diam = max( max(fwd_model.nodes(:,1:2)) -  ...
                   min(fwd_model.nodes(:,1:2)));
    gamma= 0.05*xy_diam;
end

Reg = calc_exponential_covar_prior( fwd_model, gamma);
Reg= 0.5*(Reg+Reg'); % calculation should be symmetric, but is slightly off.

Reg= inv(Reg); % FIXME - not part of inverse

try 
    noser_exp= inv_model.exponential_covar_prior.noser_exp;
catch
    noser_exp= 0;
end

if noser_exp>0
    
    J = calc_jacobian( mk_image( inv_model ));
% If J is too big, then the jacobian includes other terms (like movement) remove them
% TODO: cruft code
    n_elem = size(inv_model.fwd_model.elems,1);
    if size(J,2) > n_elem; J = J(:,1:n_elem); end
    l_prior= size(J,2);

    % Reg is spdiags(diag(J'*J),0, l_prior, l_prior);
    diag_col= sum(J.^2,1)';
    RegN = spdiags( diag_col.^noser_exp/2, 0, l_prior, l_prior);

    Reg= RegN * Reg * RegN;
end

% Calculate exponential LPF Filter
% parameter is gamma (normally 0.1)
function Reg= calc_exponential_covar_prior( fwd_model, gamma)
   [rad,ctr]= get_elem_rad_ctr( fwd_model );

   n_elem= size(ctr,1);
   oo= ones(n_elem,1);
   radh = rad/2;
   Reg=zeros(n_elem,n_elem);
   for i=1:size(ctr,1)
      ctr_i = sqrt( sum( (ctr - oo*ctr(i,:)).^2 , 2));
      Reg_i= integ_fn(-radh,radh,ctr_i-radh(i),ctr_i+radh(i), gamma);
      Reg(:,i)= Reg_i .*(Reg_i>1e-4);
   end
   Reg=sparse(Reg);

function [rad,elem_ctr]= get_elem_rad_ctr( fwd_model );
   pp= aa_fwd_parameters( fwd_model);
   if     pp.n_dims==2 % in 2d A=pi*r^2
      rad= sqrt(pp.VOLUME/pi);
   elseif pp.n_dims ==3 % in 3D V=4/3*pi*r^3
      rad= (pp.VOLUME*3/4/pi).^(1/3);
   elseif pp.n_dims ==1 % in 1D V=2*r
      rad= pp.VOLUME/2;
   else 
      error('dont know what to do with n_dims ==%d',pp.n_dims);
   end
%  rad =   rad* ones(1,pp.n_elem); % copy to matrix

   node_map = fwd_model.nodes(pp.ELEM,:);
   elem_ctr = mean( reshape( node_map, pp.n_dims+1, [], pp.n_dims), 1);
   elem_ctr = squeeze( elem_ctr);
% show_fem(fwd_model); hold on; plot(elem_ctr(:,1),elem_ctr(:,2),'*'); hold off

%  ctr= zeros(pp.n_elem);
%  for i=1:pp.n_dims
%     v= elem_ctr(:,i)*ones(1,pp.n_elem); % make square
%     ctr = ctr + (v - v').^2;
%  end
%  ctr = sqrt(ctr);

% calculate the exponential integral over a space x1,x2
%  given gamma, x1, x2, y1, y2
function fi= integ_fn(x1i,x2i,y1i,y2i, gamma)
   i_gam = 1/gamma;

   x1= min(x1i,x2i); x2= max(x1i,x2i);
   y1= min(y1i,y2i); y2= max(y1i,y2i);

   Dx= abs(x1-x2); Dy= abs(y1-y2);

   xa1= x1; 
   xa2= max(x1,min(x2,y1)); xb1= xa2;
   xb2= max(x1,min(x2,y2)); xc1= xb2;
   xc2= x2;

   RA= exp(-i_gam*y1).*(1-exp(-i_gam*Dy));
   RC= exp( i_gam*y2).*(1-exp(-i_gam*Dy));
   fi = 1./Dx./Dy/i_gam^2 .* ( ... 
            RA.*(exp( i_gam*xa2) - exp( i_gam*xa1)) - ...
            RC.*(exp(-i_gam*xc2) - exp(-i_gam*xc1)) + ...
            2*i_gam*(xb2-xb1) - ...
            exp(-i_gam*y2).*(exp( i_gam*xb2) - exp( i_gam*xb1)) + ...
            exp( i_gam*y1).*(exp(-i_gam*xb2) - exp(-i_gam*xb1))...
           );
   
   
