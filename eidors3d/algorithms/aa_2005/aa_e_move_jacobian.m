function J= aa_e_move_jacobian( fwd_model, img)
% AA_E_MOVE_JACOBIAN: J= aa_e_move_jacobian( fwd_model, img)
% Calculate Jacobian Matrix for EIT, based on conductivity
%   change and movement of electrodes
% J         = Jacobian matrix
% fwd_model = forward model
%
% fwd_model.conductivity_jacobian = fcn
%
% fwd_model.normalize_measurements if param exists, calculate
%                                  a Jacobian for normalized
%                                  difference measurements
% img = image background for jacobian calc

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: aa_e_move_jacobian.m,v 1.13 2007-08-29 09:04:03 aadler Exp $

pp= aa_fwd_parameters( fwd_model );
delta= 1e-6; % tests indicate this is a good value
             % too high and J is not linear, too low and numerical error

if isfield(fwd_model,'conductivity_jacobian')
   Jc= feval(fwd_model.conductivity_jacobian, fwd_model, img );
else
   Jc= conductivity_jacobian_perturb( pp, delta, img );
end

Jm= movement_jacobian( pp, delta, img );
J= [Jc,Jm];

% calculate normalized Jacobian
if pp.normalize
   data= fwd_solve( img );
   J= J ./ (data.meas(:)*ones(1,size(J,2)));
end


% Calculate the conductivity jacobian based on a perturbation
% of each element by a delta 
% Relative error mean(mean(abs(J-Jx)))/ mean(mean(abs(J)))
%   10^-2   0.00308129369015
%   10^-3   0.00030910899216
%   10^-4   0.00003092078190
%   10^-5   0.00000309281790
%   10^-6   0.00000035468582
%   10^-7   0.00000098672198
%   10^-8   0.00000938262464
%   10^-9   0.00009144743903

function J= conductivity_jacobian_perturb( pp, delta, img );

J = zeros( pp.n_meas, pp.n_elem );

elem_data = img.elem_data;
d0= fwd_solve( img );
for i=1:pp.n_elem
   img.elem_data   = elem_data;
   img.elem_data(i)= elem_data(i) + delta;
   di= fwd_solve( img );
   J(:,i) = (1/delta) * (di.meas - d0.meas);
end

% xy-Movement Jacobian
function J= movement_jacobian( pp, delta, img )

J = zeros( pp.n_meas, pp.n_elec*pp.n_dims );

node0= img.fwd_model.nodes;
d0= fwd_solve( img );
for d= 1:pp.n_dims
   for i= 1:pp.n_elec
      idx= img.fwd_model.electrode(i).nodes;

      img.fwd_model.nodes( idx, d)= node0(idx,d) + delta;
      di= fwd_solve( img );
      img.fwd_model.nodes( idx, d)= node0(idx,d);

      J_idx = pp.n_elec*(d-1) + i;
      J(:,J_idx) = (1/delta) * (di.meas - d0.meas);
   end
end
