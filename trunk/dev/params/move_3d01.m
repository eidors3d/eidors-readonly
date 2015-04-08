% Simulate 3D movement $Id$

noiselev = .1;
movement = 2;

% Generate eidors 3D finite element model
mdl3dim = mk_common_model( 'n3r2', [16 2]);
mdl3dim.fwd_model.nodes(:,3) = mdl3dim.fwd_model.nodes(:,3)/3;

img = mk_image(mdl3dim);
vh = fwd_solve( img );

load('datacom.mat','A','B');
img.elem_data(A) = 1.2; 
img.elem_data(B) = 0.8;
node0 = img.fwd_model.nodes;  % Node variable before movement
node1 = node0;                % Node variable after movement
z_axis = node1(:,3);

exaggeration = 10;
movement= 0.01*movement;
% Do a 3D twist - exaggerated for clearer illustration of distortion
node1(:,1) = node0(:,1).*(1 + exaggeration*movement*z_axis);
node1(:,2) = node0(:,2).*(1 + exaggeration*movement*(1-z_axis));

img.fwd_model.nodes = node1;
show_fem( img );
xlabel('x'); ylabel('y');
view(-44,22)

print_convert  move_3d01.png '-density 75'

% Do a 3D twist - we'll actually use this one
centre = 1- movement/2;
node1(:,1) = node0(:,1).*(centre + movement*z_axis);
node1(:,2) = node0(:,2).*(centre + movement*(1-z_axis));

% Solve inhomogeneous forward problem with movements and normal noise.
img.fwd_model.nodes = node1;
vi = fwd_solve( img );
noise = noiselev*std( vh.meas - vi.meas )*randn( size(vi.meas) );
vi.meas = vi.meas + noise;
move = node1 - node0;
