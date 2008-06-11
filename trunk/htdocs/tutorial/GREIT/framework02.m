% image prior $Id: framework02.m,v 1.2 2008-06-11 14:56:22 aadler Exp $
load Jacobian J

% inefficient code - but for clarity
diagJtJ = diag(J'*J);

R= spdiags( diagJtJ,0, length(diagJtJ), length(diagJtJ));

save ImagePrior R
