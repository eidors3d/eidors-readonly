% fwd_model $Id: framework02.m,v 1.1 2008-05-19 21:19:23 aadler Exp $
load Jacobian J

% inefficient code - but for clarity
diagJtJ = diag(J'*J);
