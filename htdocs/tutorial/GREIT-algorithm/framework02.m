% image prior $Id$
load Jacobian J

% inefficient code - but for clarity
diagJtJ = diag(J'*J);

R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

save ImagePrior R
