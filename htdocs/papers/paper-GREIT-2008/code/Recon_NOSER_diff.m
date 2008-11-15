function Recon_NOSER_diff( savename );
% Reconstruct GREIT images using NOSER algorithm (difference)
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM,map] = calc_NOSER_RM;
   normalize_flag = 0;
   save(savename, 'RM','normalize_flag');

function [RM,map] = calc_NOSER_RM
   [J,map] = calc_jacobian_mdl;
   RM = zeros(size(J'));
 
   % Remove space outside FEM model
   J= J(:,map);
   % inefficient code - but for clarity
   diagJtJ = diag(J'*J);
   R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

   hp = 2.0;
   RM(map,:)= (J'*J + hp^2*R)\J';
