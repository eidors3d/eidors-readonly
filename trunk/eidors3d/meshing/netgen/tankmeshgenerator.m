% Script to generate tank model interactively
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
display('EIDORS 3D tank generation script');
display('You have to answe lots of tedious questions');
display('Someone please make a GUI!');


CorR = upper(input('Circular electrodes or rectangular? [C/R]','s'));
tank_radius = input('Tank radius? ');
tank_height = input('Tank height? ');
display('Number of electrodes on each plane is 2^k? ');
elec_per_plane_index=input('Input the index k? ');
elecs_per_plane= 2^elec_per_plane_index;

display(sprintf('Thats %d electrodes per plane? ',elecs_per_plane));
no_of_planes = input('Number of planes? ');


nelec = no_of_planes*elecs_per_plane;
display(sprintf('Ok so you have %d electrodes in total? ',nelec));
if nelec < 9
  display('Thats not very many!');
elseif nelec <17
  display('Come on this is 3D EIT. Its not the 1980s you know!');
elseif nelec <65
  display('Ok sounds like enough!');
else
  display('Awesome!');
end
first_plane_starts = input('Height of centre of first plane? ');
height_between_centres = input('Height between centres? ');
if CorR=='C'
   electrode_radius= input('Electrode radius? ');
else
   electrode_height= input('Electrode height? ');
   electrode_width = input('Electrode width? ');
end

% Need some sanity checks here on the data!

fnstem = input('Enter file base file name for .geo and .vol file (ext is added):','s');

geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
[fid,mess]=fopen(geofn,'w');

if CorR =='R'

% Rectangular case 
[fid,mess]=fopen(geofn,'w');
fprintf(fid,'#Automatically generated cylinder surface mesh with rectangular electrodes\n');
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
%    fprintf(fid,'solid cyl%d = bigcyl and  cylinder %(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f);\n',kcyl,-x,-y,z,x,y,z,electrode_radius );
 end;
end;


for k= 1:nelec
  fprintf(fid,'solid cyl%d = bigcyl    and rod%d; \n',k,k);
end


fprintf(fid,'tlo bigcyl;\n');
for k=1:nelec
%  fprintf(fid,'tlo rod%d -col=[1,0,0];\n ',k);
  fprintf(fid,'tlo cyl%d cyl -col=[1,0,0];\n ',k);
end

else

% Circular case




fprintf(fid,'#Automatically generated cylinder surface mesh with circular electrodes\n');
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
    writengcylrod(fid,sprintf('rod%d',kel),[x,y,z],dirn,electrode_radius, 0.5*tank_radius);	 	 
%    fprintf(fid,'solid cyl%d = bigcyl and  cylinder %(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f);\n',kcyl,-x,-y,z,x,y,z,electrode_radius );
 end;
end;


for k= 1:nelec
  fprintf(fid,'solid cyl%d = bigcyl      and rod%d; \n',k,k);
end


fprintf(fid,'tlo bigcyl;\n');
for k=1:nelec
%  fprintf(fid,'tlo rod%d -col=[1,0,0];\n ',k);
  fprintf(fid,'tlo cyl%d cyl -col=[1,0,0];\n ',k);
end
end % of circular case

fclose(fid);
% Now call Netgen in batchmode to mesh this CSG file
disp('Calling Netgen. Please wait.....');
status= system(sprintf('ng -batchmode -geofile=%s  -meshfile=%s ',geofn,meshfn));
if status~=0
   error('Netgen call failed. Is netgen installed and on the search path?');
else
   display('Netgen seems to have meshed your tank ok and written it to file')
end
disp(['Now reading back data from file: ' meshfn])
[srf,vtx,fc,bc,simp,edg,mat_ind] = FEM_read_mesh_2(meshfn);
disp([meshfn ' contains ' num2str(max(fc)) ' faces'])

disp('Now I need some help finding which faces are electrodes')

disp('Do not choose and ground planes unless you have a vessel with a partly metal wall!')


% Plot wire frame equivalent of mesh model
%tetramesh(simp,vtx,'FaceColor','none','EdgeColor','cyan')%it doesn't work in matlab 5.3
set(gcf,'Name','Wire Mesh Model')
view(45,10)
hold on
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'EdgeColor','blue')
title('Surface of body: blue, Volume of body: cyan')
hidden off
axis equal image; % Tightly fit square axes around plot
mshaxs = axis; % Save present axes for use with faces
pause(3)

% Select the electrodes and ground planes (last electrode)
[elecgnd,sels,sgnd] = FEM_define_electrodes_2(srf,vtx,bc,mshaxs);




