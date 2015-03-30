% Solve 2D and 3D model $Id$

% Original target
subplot(141)
show_fem(img); view(-62,28)

% Create inverse Model: Classic
imdl= select_imdl(f_mdl, {'Basic GN dif'});
imdl.hyperparameter.value = .1;

% Classic (inverse crime) solver
img1= inv_solve(imdl, vh, vi);
subplot(142)
show_fem(img1); view(-62,28)
