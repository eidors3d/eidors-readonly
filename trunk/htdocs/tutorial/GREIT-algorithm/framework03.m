% reconst $Id$
load ImagePrior R
load Jacobian J

hp = .01;

RM= (J'*J + hp^2*R)\J';
