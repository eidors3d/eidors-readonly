function [img,map]= GREIT_NOSER_ndiff( ref_meas, reconst_meas )
% Reconstruct GREIT images using NOSER algorithm
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM,map] = calc_NOSER_RM;

   % Expand ref_meas to the full size of reconst_meas
   num_meas = size(reconst_meas,2);
   ref_meas = ref_meas * ones(1,num_meas);
   dv = ( reconst_meas - ref_meas ) ./ ref_meas; % CHANGE IS HERE:

   % reconst image
   ds = RM*dv;

   img= reshape(ds, 32,32,num_meas);

function [RM,map] = calc_NOSER_RM
   [J,map,vbkgnd] = GREIT_Jacobian_cyl;
   J = J ./ (vbkgnd*ones(1,size(J,2))); % Normalized Jacobian
   RM = zeros(size(J'));
 
   % Remove space outside FEM model
   J= J(:,map);
   % inefficient code - but for clarity
   diagJtJ = diag(J'*J);
   R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

   hp = 2.0;
   RM(map,:)= (J'*J + hp^2*R)\J';
