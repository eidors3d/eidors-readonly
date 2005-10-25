function Reg= aa_e_move_image_prior( inv_model );
% AA_E_MOVE_IMAGE_PRIOR calculate image prior
% Reg= aa_e_move_image_prior( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
%
% For image portion, we use a Laplace prior, as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself
%
% For the movmenent portion, we define a smoothness
% constraint, such that Rij = -1 for adjacent electrodes

% $Id: aa_e_move_image_prior.m,v 1.1 2005-10-25 17:18:27 aadler Exp $

% relative strengths of conductivity and movement priors
hp_cond= 1;
hp_move= 1;

pp= aa_fwd_parameters( inv_model.fwd_model );

% calc conductivity portion
RegC = laplace_image_prior( inv_model );

% calc movement portion
RegM = movement_image_prior( pp.n_dims, pp.n_elec );

% zero blocks
RegCM= sparse( pp.n_elem, pp.n_dims*pp.n_elec );

Reg= [hp_cond* RegC,           RegCM;
               RegCM', hp_move*RegM ];

% For the movmenent portion, we define a smoothness
% constraint, such that Rij = -1 for adjacent electrodes
function RegM = movement_image_prior( dims, elecs );

   % movement constraint for 1 dimention
   idx =(1:elecs)';
   im1= rem(idx-1+elecs,elecs);
   ip1= rem(idx+1,elecs); 
   mv= sparse([idx,idx,idx]+1,[im1,idx,ip1]+1,ones(elecs,1)*[-1,2,-1]);

   RegM= spalloc(dims*elecs,dims*elecs, 3*dims*elecs);

   for i=0:dims-1;
     RegM( idx+i*elecs, idx+i*elecs) = mv;
   end

