% reconst $Id$
load ImagePrior_diag_JtJ R
load GREIT_Jacobian_ng_mdl_fine J map

RM = zeros(size(J'));
J = J(:,map);

hp = .3;

RM(map,:)= (J'*J + hp^2*R)\J';

save RM_framework_example RM map
