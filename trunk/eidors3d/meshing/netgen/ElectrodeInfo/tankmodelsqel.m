% Script to make model of tank with square electrodes
% dimensions in mm
% height 180mm
% diam 140
%diam elec 9mm
tank_radius = 20
tank_height = 35
elec_per_plane_index=3
elecs_per_plane= 2^elec_per_plane_index
no_of_planes = 2
first_plane_starts = 10
height_between_centres = 15
electrode_height= 8
electrode_width =5
%tank_radius = 20
%tank_height = 35
%no_of_planes = 2
%elecs_per_plane= 8
%first_plane_starts = 10
%height_between_centres = 15
%electrode_radius= 2
eps = tank_radius*0.1;
geofn='C:\EIDORS3D\netgen\tanksqel.geo';
meshfn='tanksqel.msh';
[fid,mess]=fopen(geofn,'w');
fprintf(fid,'#Automatically generated cylinder surface mesh with square electrodes\n');
fprintf(fid,'algebraic3d\n');
fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,%6.2f;%6.2f); \n',tank_height, ...
	tank_radius);
fprintf(fid, 'solid bigcyl= plane(0,0,0;0,0,-1)\n  and  plane(0,0,%6.2f;0,0,1)\n   and  cyl;\n',tank_height);  
kel = 0;
for l=1:no_of_planes
 z = first_plane_starts + (l-1)*height_between_centres;
 for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/elecs_per_plane
    kel=kel+1;
    [x,y]=pol2cart(th,tank_radius);
    dirn = [x,y,0];
    dirn = dirn ./ norm(dirn);
    writengcuboid(fid,sprintf('rod%d',kel),[x,y,z],dirn,electrode_height,electrode_width,0.5*tank_radius);	 	 
    %fprintf(fid,'solid cyl%d = bigcyl and  cylinder %(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f);\n',kcyl,-x,-y,z,x,y,z,electrode_radius );
 end;
end;



nelec = no_of_planes*elecs_per_plane;
for k= 1:nelec
  fprintf(fid,'solid cyl%d = bigcyl \n           and rod%d;',k,k);
end


fprintf(fid,'tlo bigcyl;\n');
for k=1:nelec
%  fprintf(fid,'tlo rod%d -col=[1,0,0];\n ',k);
  fprintf(fid,'tlo cyl%d cyl -col=[1,0,0];\n ',k);
end
fclose(fid);
