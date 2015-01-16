function Reg= prior_movement( inv_model );
% PRIOR_MOVEMENT calculate image prior
% Reg= prior_movement( inv_model )
% Reg        => output regularization term
% inv_model  => inverse model struct
% Parameters:
%   inv_model.prior_movement.parameters(1) -> relative weighting
%     of movement vs image fraction of hyperparameter
%     => Default = 10
%   inv_model.prior_movement.RegC.func = Cond Reg fcn
%   inv_model.prior_movement.RegM.func = Move Reg fcn
%   either @laplace_movement_image_prior OR @tikhonov_movement_image_prior
%
% For image portion, we use a Laplace prior, as 
% -1 for each adjacent element, and 3 (in 2D) or 4 (in 3D)
% for the element itself
%
% For the movement portion, we define a smoothness
% constraint, such that Rij = -1 for adjacent electrodes
%
% If used with a dual model (ie coarse2fine mapping), ensure
%    imdl.prior_use_fwd_not_rec = 1;

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id$

% relative strengths of conductivity and movement priors
if isfield( inv_model,'prior_movement')
   hp_move= inv_model.prior_movement.parameters(1);
else
   hp_move= 10;
end

pp= fwd_model_parameters( inv_model.fwd_model );
n_elec = pp.n_elec;

% calc conductivity portion
try
   RegCfcn = inv_model.prior_movement.RegC.func;
catch
   RegCfcn = @prior_laplace;
end
try; inv_model = rmfield(inv_model, 'R_prior'); end
try; inv_model = rmfield(inv_model, 'prior_use_fwd_not_rec'); end
inv_model.RtR_prior = RegCfcn;

pp= fwd_model_parameters( inv_model.fwd_model );


RegC= calc_RtR_prior( inv_model); 
szC = size(RegC,1);

% calc movement portion
try
   fname = func2str(inv_model.prior_movement.RegM.func);
   if(strcmp(fname,'tikhonov_movement_image_prior'))
       RegM = tikhonov_movement_image_prior(pp.n_dims,pp.n_elec);
   elseif(strcmp(fname,'laplace_movement_image_prior'))
       RegM=laplace_movement_image_prior(pp.n_dims,pp.n_elec);
   end
catch
   RegM = laplace_movement_image_prior( pp.n_dims, pp.n_elec );
end
% zero blocks
RegCM= sparse( szC, pp.n_dims*pp.n_elec );

Reg= [RegC,           RegCM;
      RegCM', hp_move^2*RegM ]; 

% For the movmenent portion, we define a smoothness
% constraint, such that Rij = -1 for adjacent electrodes
%
% TODO: this approach assumes that electrodes close in number
%   are close to each other. This isn't necessarily right.
function RegM = laplace_movement_image_prior( dims, elecs );

   % movement constraint in each dimention
   idx =(0:elecs-1)';
   im1= rem(idx-1+elecs,elecs);
   ip1= rem(idx+1,elecs); 
   mv= sparse([idx,idx,idx]+1,[im1,idx,ip1]+1,ones(elecs,1)*[-1,2.1,-1]);

   RegM= spalloc(dims*elecs,dims*elecs, 3*dims*elecs);

   for i=0:dims-1;
     m_idx= idx + i*elecs + 1;
     RegM( m_idx, m_idx ) = mv;
   end
   
function RegM = tikhonov_movement_image_prior( dims, elecs);
   RegM=eye(dims*elecs);   
