% image prior $Id$
load GREIT_Jacobian_ng_mdl_fine J

% inefficient code - but for clarity
diagJtJ = diag(J'*J);

R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

save ImagePrior_diag_JtJ R
