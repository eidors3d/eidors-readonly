function Recon_NOSER_ndiff( savename );
% Reconstruct GREIT images using NOSER algorithm
%
% (C) 2008 Andy Adler. Licenced under GPL v2 or v3
% $Id$

   [RM,map] = calc_NOSER_RM;
   normalize_flag = 1;
   save(savename, 'RM','normalize_flag');

function [RM,map] = calc_NOSER_RM
   [J,map,vbkgnd] = calc_jacobian_mdl;
   J = J ./ (vbkgnd*ones(1,size(J,2))); % Normalized Jacobian
   RM = zeros(size(J'));
 
   % Remove space outside FEM model
   J= J(:,map);
   % inefficient code - but for clarity
   diagJtJ = diag(J'*J);
   R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

   hp = 2.0;
   RM(map,:)= (J'*J + hp^2*R)\J';
