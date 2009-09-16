function [mdl, mat_indices] = create_2d2l_mesh_gmsh(basename, in_rad, ext_rad, no_elec)
% Create a 2D Circular FEM using GMSH
% mdl= CREATE_GMSH_2D_CIRCLE(rad, n_elec)
%
% mdl - EIDORS forward model
% rad - model radius

% (C) 2006 Andy Adler. License: GPL version 2 or version 3
% $Id: create_gmsh_2d_circle.m 1535 2008-07-26 15:36:27Z aadler $


geo_filename = sprintf('%s.geo', basename);
geo_fid= fopen(geo_filename,'w');

% Electrodes points on external boundary
theta= linspace(0,2*pi, no_elec+1); theta(end)=[];
point_no = 1;
for th=theta
    x=ext_rad*sin(th);
    y=ext_rad*cos(th);
    z=0;
    fprintf(geo_fid,'Point(%d) = {%f,%f,%f,%f};\n',point_no, x,y,z, ext_rad/100);
    point_no = point_no + 1;
end

% Points to define internal layer
center_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f,%f,%f,%f};\n',center_no, 0, 0, 0, 1);
point_no = point_no + 1;
inpoint1_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f,%f,%f,%f};\n',inpoint1_no, in_rad, 0,0, ext_rad/30);
point_no = point_no + 1;
inpoint2_no = point_no;
fprintf(geo_fid,'Point(%d) = {%f,%f,%f,%f};\n',inpoint2_no, -in_rad, 0,0, ext_rad/30);
point_no = point_no + 1;

% External lines / electrodes
line_no = 1;
for i = 1:no_elec
    start_point = i;
    end_point= i+1;
    if end_point > no_elec
        end_point = 1;
    end
    fprintf(geo_fid,'Line(%d) = {%d, %d};\n', line_no, start_point, end_point);
    line_no = line_no + 1;
end

extloop_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Line Loop(%d) = {', extloop_no);
for i = 1:(no_elec-1)
    fprintf(geo_fid,'%d,', i);
end
fprintf(geo_fid, '%d};\n', no_elec);

% Internal circle and loop
circle1_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Circle(%d) = {%d, %d, %d};\n', circle1_no, inpoint1_no, ...
    center_no, inpoint2_no);
circle2_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Circle(%d) = {%d, %d, %d};\n', circle2_no, inpoint2_no, ...
    center_no, inpoint1_no);

inloop_no = line_no;
line_no = line_no + 1;
fprintf(geo_fid,'Line Loop(%d) = {%d,%d};\n', inloop_no, circle1_no, ...
    circle2_no);

fprintf(geo_fid, 'Plane Surface(%d) = {%d, %d};\n', line_no, extloop_no, ...
    inloop_no);
line_no = line_no + 1;

fprintf(geo_fid, 'Plane Surface(%d) = {%d};\n', line_no, inloop_no);
line_no = line_no + 1;

fclose(geo_fid);

% Call Gmsh 
disp('Calling Gmsh. Please wait ...');
call_gmsh( geo_filename);

msh_filename = sprintf('%s.msh', basename);

disp(['Now reading back data from file: ' msh_filename])
[mdl, mat_indices]= gmsh_mk_fwd_model( msh_filename, ...
    'Gmsh based 2 layer circural model', ones(no_elec,3) );


