% Script to make model of tank in E&EE
% dimensions in mm
% height 180mm
% diam 140
%diam elec 9mm
tank_radius = 90
tank_height = 140
elecs_per_plane= 16
no_of_planes = 5
first_plane_starts = 30
height_between_centres = 30
electrode_radius= 9
%tank_radius = 20
%tank_height = 35
%no_of_planes = 2
%elecs_per_plane= 8
%first_plane_starts = 10
%height_between_centres = 15
%electrode_radius= 2
eps = 1e-8;
geofn='tank.geo';
meshfn='tank.msh';
[fid,mess]=fopen(geofn,'w');
fprintf(fid,'#Automatically generated cylinder surface mesh with round electrodes\n');
fprintf(fid,'algebraic3d\n');
fprintf(fid,'solid cy=cylinder(0,0,0;0,0,%6.2f;%6.2f) \n',tank_height, ...
	tank_radius);
fprintf(fid,' and plane(0,0,0;0,0,-1)\n  and  plane(0,0,%6.2f;0,0,1);\n',tank_height);  
kelec = 0;
for l=1:no_of_planes
 z = first_plane_starts + (l-1)*height_between_centres;
 for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/ elecs_per_plane
    kelec=kelec+1;
    [x,y]=pol2cart(th,tank_radius);
	 	 
    fprintf(fid,'solid el%d = sphere(%6.2f,%6.2f,%6.2f;%6.2f);\n',kelec,x+eps,y+eps,z+eps,electrode_radius );
 end;
end;
nelec = no_of_planes*elecs_per_plane;
fprintf(fid,'solid el = ');
for k=1:(nelec-1)
  fprintf(fid,'el%d or ',k);
end
fprintf(fid,'el%d;\n',nelec);
fprintf(fid,'#cylinder with spheres scooped out\n');
fprintf(fid,'solid body = cy and not el;\n');
fprintf(fid,'#the bits which were scooped out\n');
fprintf(fid,'solid e1 = cy and el;\n');
fprintf(fid,'#Making two top level objects means they both get meshed.\n');
fprintf(fid,'tlo body -transparent;\n');
fprintf(fid,'tlo e1;\n');     
fclose(fid);
% Now call Netgen in batchmode to mesh this CSG file
disp('Calling Netgen. Please wait.....');
status= system(sprintf('ng -batchmode -geofile=%s  -meshfile=%s ',geofn,meshfn));
if status~=0
   error('Netgen call failed. Is netgen installed and on the search path?');
end
% Check of we can read the mesh in
[vtx,simp,surf] = readngvol(meshfn);
size(vtx)
size(simp)
size(surf)