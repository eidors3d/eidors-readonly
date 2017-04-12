% Solve 2.5D
% Create inverse Model: Classic
imdl.hyperparameter.value = .1;
imdl.reconst_type = 'difference';
imdl.type = 'inv_model';

% Classic and 2.5D (inverse crime) solver
vh.type = 'data';
vi.type = 'data';
img2 = inv_solve(imdl, vh, vi);
imdl.fwd_model.jacobian = @jacobian_adjoint_2p5d_1st_order;
img25 = inv_solve(imdl, vh, vi);

clf;
subplot(131); show_fem(fimg); title('model'); axis off; axis square;
subplot(132); show_fem(img2); title('2D'); axis off; axis square;
subplot(133); show_fem(img25);title('2.5D'); axis off; axis square;
print_convert two_and_half_d03a.png
