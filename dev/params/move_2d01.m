% Generate simulation data without noise and standard reconstruction
% Create circular FEM - creates a eidors_mdl type inv_model.
mdlc = mk_common_model('c2c');

f_img = mk_image( mdlc, 1);
vh = fwd_solve( f_img );

% Hard coded values here represent local inhomogeneities
f_img.elem_data([75,93,94,113,114,136]) = 1.2;
f_img.elem_data([105,125,126,149,150,174]) = 0.8;

% Simulate node movements - shrink x, stretch y by 1% of model diameter
% node0 before, node1 after movement

movement = [1-0.01 0; 0 1+0.01];
node0 = f_img.fwd_model.nodes;
node1 = node0*movement;
f_img.fwd_model.nodes = node1;

% Solve inhomogeneous forward problem with movements and normal noise
% 1% of standard deviation of signal
vi = fwd_solve( f_img );
move = node1 - node0;

% Plot FEM with conductivities and movement vectors.
show_fem_move( f_img, move, 20 );
print_convert move_2d01.png '-density 75'
