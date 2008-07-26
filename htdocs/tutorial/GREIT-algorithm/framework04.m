% reconst $Id$
load ImagePrior_diag_JtJ R
load GREIT_Jacobian_ng_mdl_fine J

hp = .01;

RM= (J'*J + hp^2*R)\J';
