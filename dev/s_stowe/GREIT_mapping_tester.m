body_geometry.cone = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) 0.6];
    electrode_geometry{i}.sphere.radius = 0.05;
end
for i = 17:32
    electrode_geometry{i}.sphere.center = [0.75*cos(theta(i-16)) 0.75*sin(theta(i-16)) 0.4];
    electrode_geometry{i}.sphere.radius = 0.05;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
figure(1)
show_fem(fmdl,[0,1])
skip4 = {32,1,[0,5],[0,5],{'no_meas_current_next1'},1};
[fmdl.stimulation,fmdl.meas_select] = mk_stim_patterns(skip4{:});
% Set GREIT parameters
% vopt
vopt.imgsz = [32 32];
vopt.zvec = linspace(0.1,0.9,9);
vopt.square_pixels = true;
opt.target_plane = 0.25;

% GREIT 3D with 2x16 electrode belt
[imdl,opt.distr] = GREIT3D_distribution(fmdl, vopt);

imdl3= mk_GREIT_model(imdl, 0.5, [20], opt);

% Issue: imdl3 give a RM of all NAN
