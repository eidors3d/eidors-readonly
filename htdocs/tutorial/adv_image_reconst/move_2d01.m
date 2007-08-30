% Set eidors default colours
calc_colours('backgnd',[.9,.9,.9]);

% Generate simulation data without noise and standard reconstruction
% Create circular FEM - creates a eidors_mdl type inv_model.
mdlc = mk_common_model('c2c');

% Instantiate a homogeneous forward model.
% ref_level is 1 since we use ones( ).
sigma = ones( size(mdlc.fwd_model.elems,1) ,1);
% Create the eidors_obj of type image.
f_img = eidors_obj('image', 'homogeneous image', ...
    'elem_data', sigma, ...
    'fwd_model', mdlc.fwd_model,...
    'ref_level', 1);
% Solve homogeneous forward problem
vh = fwd_solve( f_img );

% Hard coded values here represent local inhomogeneities
sigma([75,93,94,113,114,136]) = 1.2;
sigma([105,125,126,149,150,174]) = 0.8;
f_img.elem_data = sigma;

% Simulate node movements - shrink x, stretch y by 1% of model diameter
% node0 before, node1 after movement
movement = [1-0.01 0; 0 1+0.01];
node0 = f_img.fwd_model.nodes;
node1 = node0*movement;
f_img.fwd_model.nodes = node1;

% Solve inhomogeneous forward problem with movements and normal noise
% 1% of standard deviation of signal
vi = fwd_solve( f_img );
noise = 0*std( vh.meas - vi.meas )*randn( size(vi.meas) );
vi.meas = vi.meas + noise;
move = node1 - node0;

% Plot FEM with conductivities and movement vectors.
show_fem_move( f_img, move, 20 );
print -r75 -dpng move_2d01.png
