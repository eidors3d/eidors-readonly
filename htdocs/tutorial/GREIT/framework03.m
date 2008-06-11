% reconst $Id: framework03.m,v 1.1 2008-06-11 14:56:22 aadler Exp $
load ImagePrior R
load Jacobian J

hp = .01;

RM= (J'*J + hp^2*R)\J';
