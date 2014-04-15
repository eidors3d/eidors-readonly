

%CASE  1 %%%%
disp('#### 01 ####');clear;
body_geometry.cylinder = struct;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models01.png


%CASE  2 %%%%
disp('#### 02 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
    electrode_geometry{i}.sphere.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models02.png


%CASE  3 %%%%
disp('#### 03 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models03.png


%CASE  4 %%%%
disp('#### 04 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
    electrode_geometry{i}.keep_material_flag = 1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models04.png


%CASE  5 %%%%
disp('#### 05 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
    electrode_geometry{i}.enter_body_flag = 1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models05.png


%CASE  6 %%%%
disp('#### 06 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
    electrode_geometry{i}.keep_material_flag = 1;
    electrode_geometry{i}.enter_body_flag = 1;                
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models06.png


%CASE  7 %%%%
disp('#### 07 ####');clear;
body_geometry.cylinder = struct;
body_geometry.sphere.center = [0 0 1];
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models07.png


%CASE  8 %%%%
disp('#### 08 ####');clear;
body_geometry.cylinder  = struct;
body_geometry.sphere(1) = struct;  
body_geometry.sphere(2).center = [0 0 1];         
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);   
show_fem( fmdl );
print_convert netgen_geometric_models08.png


%CASE  9 %%%%
disp('#### 09 ####');clear;
body_geometry.intersection.cylinder(1) = struct;
body_geometry.intersection.cylinder(2).radius     = 0.5;
body_geometry.intersection.cylinder(2).complement_flag = 1;   
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models09.png


%CASE 10 %%%%
disp('#### 10 ####');clear;
body_geometry.intersection(1).sphere(1).radius     = 0.5;
body_geometry.intersection(1).sphere(1).center     = [0 0 2];
body_geometry.intersection(1).sphere(1).complement_flag = 1;
body_geometry.intersection(1).sphere(2).center     = [0 0 2];
body_geometry.intersection(2).cylinder(1).top_center = [0 0 2];
body_geometry.intersection(2).cylinder(2).radius     = 0.5;
body_geometry.intersection(2).cylinder(2).top_center = [0 0 2];
body_geometry.intersection(2).cylinder(2).complement_flag = 1;   
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models10.png


%CASE 11 %%%%
disp('#### 11 ####');clear;
body_geometry.intersection.union(1).sphere.radius = 0.5;
body_geometry.intersection.union(1).sphere.center = [0 0 2];
body_geometry.intersection.union(1).cylinder.radius = 0.5;
body_geometry.intersection.union(1).cylinder.top_center = [0 0 2];
body_geometry.intersection.union(1).complement_flag = 1;
body_geometry.intersection.union(2).sphere.center = [0 0 2];
body_geometry.intersection.union(2).cylinder.top_center = [0 0 2]; 
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models11.png


%CASE 12 %%%%
disp('#### 12 ####');clear;
body_geometry.cone = struct; 
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [0.85*cos(theta(i)) 0.85*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.65*cos(theta(i)) 0.65*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models12.png


%CASE 13 %%%%
disp('#### 13 ####');clear;
body_geometry.cone = struct; 
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) 0.5];
    electrode_geometry{i}.sphere.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models13.png


%CASE 14 %%%%
disp('#### 14 ####');clear;
body_geometry.cone(1).top_center = [0 0 1.5];
body_geometry.cone(1).bottom_center = [0 0 0.5];
body_geometry.cone(2).top_center = [0 0 -1.5];
body_geometry.cone(2).bottom_center = [0 0 -0.5];
body_geometry.cylinder.top_center    = [0, 0, 0.5];
body_geometry.cylinder.bottom_center = [0, 0, -0.5];
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) 1.0];
    electrode_geometry{i}.sphere.radius = 0.1;
    electrode_geometry{i + n_elect}.sphere.center = [cos(theta(i)) sin(theta(i)) 0];
    electrode_geometry{i + n_elect}.sphere.radius = 0.15;
    electrode_geometry{i + 2*n_elect}.sphere.center = [0.75*cos(theta(i)) 0.75*sin(theta(i)) -1.0];
    electrode_geometry{i + 2*n_elect}.sphere.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models14.png


%CASE 15 %%%%
disp('#### 15 ####');clear;
body_geometry.ortho_brick = struct;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models15.png


%CASE 16 %%%%
disp('#### 16 ####');clear;
body_geometry.intersection.ortho_brick.opposite_corner_a = [0 0 0];
body_geometry.intersection.ortho_brick.opposite_corner_b = [5 5 4];
for i = 1:4; 
    for j = 1:4; 
        body_geometry.intersection.cylinder(i,j).radius = 0.15;
        body_geometry.intersection.cylinder(i,j).top_center = [i, j, 4];
        body_geometry.intersection.cylinder(i,j).bottom_center = [i, j, 2];
        body_geometry.intersection.cylinder(i,j).complement_flag = 1;
    end; 
end;
fmdl = ng_mk_geometric_models(body_geometry);    
show_fem( fmdl );
print_convert netgen_geometric_models16.png


%CASE 17 %%%%
disp('#### 17 ####');clear;
body_geometry.intersection.ortho_brick.opposite_corner_a = [0 0 0];
body_geometry.intersection.ortho_brick.opposite_corner_a = [0 0 0];
body_geometry.intersection.ortho_brick.opposite_corner_b = [5 5 4];
for i = 1:4; 
    for j = 1:4; 
        body_geometry.intersection.cylinder(i, j).radius = 0.15;
        body_geometry.intersection.cylinder(i, j).top_center    = [i, j, 4];
        body_geometry.intersection.cylinder(i, j).bottom_center = [i, j, 2];
        body_geometry.intersection.cylinder(i, j).complement_flag = 1;
        electrode_geometry{i, j, 1}.cylinder.radius        = 0.2;
        electrode_geometry{i, j, 1}.cylinder.top_center    = [i, j, 3.1];
        electrode_geometry{i, j, 1}.cylinder.bottom_center = [i, j, 2.9];
        electrode_geometry{i, j, 2}.cylinder.radius        = 0.2;
        electrode_geometry{i, j, 2}.cylinder.top_center    = [i, j, 2.2];
        electrode_geometry{i, j, 2}.cylinder.bottom_center = [i, j, 2.0];
    end; 
end;
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models17.png


%CASE 18 %%%%
disp('#### 18 ####');clear;
body_geometry.parallelepiped  = struct;
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models18.png


%CASE 19 %%%%
disp('#### 19 ####');clear;
body_geometry.parallelepiped.vertex   = [ 0;  0;  0];
body_geometry.parallelepiped.vector_a = [ 1;  1;  0];
body_geometry.parallelepiped.vector_b = [ 0;  1;  1];
body_geometry.parallelepiped.vector_c = [ 1;  0;  1];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models19.png


%CASE 20 %%%%
disp('#### 20 ####');clear;
body_geometry.intersection.ortho_brick.opposite_corner_a = [-15, -15, 0];
body_geometry.intersection.ortho_brick.opposite_corner_b = [15, 15, 5];
body_geometry.intersection.half_space.point = [0, 0, 5];
body_geometry.intersection.half_space.outward_normal_vector = [-1, -1, 5];

fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models20.png


%CASE 21 %%%%
disp('#### 21 ####');clear;
body_geometry.ellipsoid.axis_a = [1 0 0];
body_geometry.ellipsoid.axis_b = [0 2 0];
body_geometry.ellipsoid.axis_c = [0 0 3];
fmdl = ng_mk_geometric_models(body_geometry);   
show_fem( fmdl );
print_convert netgen_geometric_models21.png


%CASE 22 %%%%
disp('#### 22 ####');clear;
body_geometry.ellipsoid.axis_a = [1 0 0];
body_geometry.ellipsoid.axis_b = [0 1 1];
body_geometry.ellipsoid.axis_c = [0 -2 2];
fmdl = ng_mk_geometric_models(body_geometry);   
show_fem( fmdl );
print_convert netgen_geometric_models22.png


%CASE 23 %%%%
disp('#### 23 ####');clear;
body_geometry.elliptic_cylinder.top_center = [0, 0, 10];
body_geometry.elliptic_cylinder.bottom_center = [0, 0, 0];           
body_geometry.elliptic_cylinder.axis_a = [1 0 0];
body_geometry.elliptic_cylinder.axis_b = [0 2 0];  
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models23.png


%CASE 24 %%%%
disp('#### 24 ####');clear;
body_geometry.elliptic_cylinder.top_center = [0, 5, 5];
body_geometry.elliptic_cylinder.bottom_center = [0, 0, 0];           
body_geometry.elliptic_cylinder.axis_a = [1 0 0];
body_geometry.elliptic_cylinder.axis_b = [0 -2 2];  
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models24.png


%CASE 25 %%%%
disp('#### 25 ####');clear;
body_geometry.body_of_revolution = struct;
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models25.png


%CASE 26 %%%%
disp('#### 26 ####');clear;
body_geometry.body_of_revolution.points   = [1 1; 1 2; 2 1.5; 2 1];
body_geometry.body_of_revolution.segments = [1 2; 2 3; 3 4; 4 1];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models26.png


%CASE 27 %%%%
disp('#### 27 ####');clear;
n_points = 24;
theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
body_geometry.body_of_revolution.points   = 2 + [sin(theta) cos(theta)];
body_geometry.body_of_revolution.segments = [(1:n_points)' [(2:n_points) 1]'];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models27.png


%CASE 28 %%%%
disp('#### 28 ####');clear;
n_points = 24;
theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
body_geometry.body_of_revolution.points   = 2 + [sin(theta) cos(theta)];
body_geometry.body_of_revolution.segments = [(1:2:n_points)' (2:2:n_points)' [(3:2:n_points) 1]'];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models28.png


%CASE 29 %%%%
disp('#### 29 ####');clear;
body_geometry{1}.cylinder(1).radius        = 0.5;
body_geometry{1}.cylinder(1).top_center    = [0 0 0.75];
body_geometry{1}.cylinder(1).bottom_center = [0 0 0.25];
body_geometry{1}.name                      = 'Object';           
body_geometry{2}.cylinder(2).radius        = 1;
body_geometry{2}.name                      = 'Tank';
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.cylinder.top_center    = [1.03*cos(theta(i)) 1.03*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.bottom_center = [0.97*cos(theta(i)) 0.97*sin(theta(i)) 0.5];
    electrode_geometry{i}.cylinder.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models29.png


%CASE 30 %%%%
disp('#### 30 ####');clear;
body_geometry{1}.sphere.radius     = 0.25;
body_geometry{1}.sphere.center     = [0 0 0.5];
body_geometry{1}.name              = 'Sphere';
body_geometry{2}.cylinder.radius   = 1;
body_geometry{2}.name              = 'Tank';           
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
    electrode_geometry{i}.sphere.radius = 0.1;
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models30.png


%CASE 31 %%%%
disp('#### 31 ####');clear;
n_sphere = 8;
theta = linspace(0, 2*pi, n_sphere+1); theta(end) = [];   
for i = 1:n_sphere
    body_geometry{i}.sphere.radius   = 0.2;
    body_geometry{i}.sphere.center   = [0.65*cos(theta(i)) 0.65*sin(theta(i)) 0.5];  
    body_geometry{i}.max_edge_length = 0.025*(1 + rem(i,2));
    body_geometry{i}.name            = sprintf('Sphere%d', i);  
end        
body_geometry{n_sphere+1}.cylinder.radius = 1;
body_geometry{n_sphere+1}.name            = 'Tank';  
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
    electrode_geometry{i}.sphere.radius = 0.1;
    electrode_geometry{i}.max_edge_length = 0.025*(1 + rem(i,2));
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models31.png


%CASE 32 %%%%
disp('#### 32 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    electrode_geometry{i}.point = [cos(theta(i)) sin(theta(i)) 0.5];
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models32.png


%CASE 33 %%%%
disp('#### 33 ####');clear;
body_geometry.cylinder = struct;
n_elect = 16;
theta = linspace(0, 2*pi, n_elect+1); theta(end) = [];
for i = 1:n_elect
    if (rem(i,2))
        electrode_geometry{i}.point = [cos(theta(i)) sin(theta(i)) 0.5];
        electrode_geometry{i}.name  = sprintf('Point_Electrode%d', ceil(i/2));
    else
        electrode_geometry{i}.sphere.center = [cos(theta(i)) sin(theta(i)) 0.5];
        electrode_geometry{i}.sphere.radius = 0.1;
        electrode_geometry{i}.name          = sprintf('Circular_Electrode%d', floor(i/2));
    end
end
fmdl = ng_mk_geometric_models(body_geometry, electrode_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models33.png


%CASE 34 %%%%
disp('#### 34 ####');clear;
body_geometry.body_of_extrusion = struct;
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models34.png


%CASE 35 %%%%
disp('#### 35 ####');clear;
body_geometry.body_of_extrusion.path_points   = [0 0 0; 0.25 0 1; 0.25 0 2; 0.25 0 3; 0 0 4];
body_geometry.body_of_extrusion.path_segments = [1 2; 2 3; 3 4; 4 5];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry);
show_fem( fmdl );
print_convert netgen_geometric_models35.png


%CASE 36 %%%%
disp('#### 36 ####');clear;
n_points = 16;
theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];
body_geometry.body_of_extrusion.profile_points   = 0.2*(2 + [0.75*sin(theta) cos(theta)]);
body_geometry.body_of_extrusion.profile_segments = [(1:n_points)' [(2:n_points) 1]'];
n_points = 32;
theta = linspace(0, 2*pi, n_points+1)'; theta(end) = [];          
body_geometry.body_of_extrusion.path_points   = 1*(2 + [sin(theta) 1.5*cos(theta) zeros(n_points, 1)]);
body_geometry.body_of_extrusion.path_segments = [(1:n_points)' [(2:n_points) 1]'];
body_geometry.body_of_extrusion.vector_d      = [0; 0; 1];
body_geometry.max_edge_length = 0.15;
fmdl = ng_mk_geometric_models(body_geometry); 
show_fem( fmdl );
print_convert netgen_geometric_models36.png
