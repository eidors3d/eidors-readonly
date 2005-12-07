function [tank_mdl,centres] = create_tank_mesh_ng( tank_radius, tank_height,CorR, log2_electrodes_per_plane, no_of_planes,first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem)
%function tank_mod = create_tank_mesh_ng( tank_radius, tank_height, circ_or_rec, log2_electrodes_per_plane, no_of_planes,..
%number_of_planes, first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem)
%
% Parameters  tank_radius, tank_height, 
%   CorR= 'C' for circular 'R' for rectangular
%  log2_electrodes_per_plane  that is 2^this electrodes per plane
%, no_of_planes,..
%  number_of_planes, 
%  first_plane_starts  z coordinate were centre of first plane starts
%, height_between_centres   of electrode planes , 
% electrode_width, electrode_height , the width is just radius if 'R'
%  fnstem the file name used for saving netgen files
% Function to generate tank model
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
% Revised for new version 3.0 structure WRBL 05/12/2005
%Made in to function WRBL 6/5/2005

elecs_per_plane= 2^log2_electrodes_per_plane;


nelec = no_of_planes*elecs_per_plane;
if CorR=='C'
   electrode_radius=electrode_width; 
end


% Need some sanity checks here on the data!


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
    centres(kel,:)= [x,y,z]; % keep the centres
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
    centres(kel,:)= [x,y,z]; % keep the centres
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
while( 1 )
   disp('Calling Netgen. Please wait.....');
   status= system(sprintf( ...
        'ng -batchmode -geofile=%s  -meshfile=%s ',geofn,meshfn));
   if status==0; break; end

   fprintf([ ...
    'Netgen call failed. Is netgen installed and on the search path?\n' ...
    'If you are running under windows, I can attempt to create\n' ...
    'a batch file to access netgen.\n' ...
    'Please enter the directory in which to find netgen.' ...
    'If you don''t have a copy, press Ctrl-C to break, and' ...
    'see http://www.hpfem.jku.at/netgen/ for download\n\n' ...
    ]);
   netgen_path = input('netgen_path? ','s');
   fid= fopen('ng.bat','w');
   fprintf(fid,'set TCL_LIBRARY=%s/lib/tcl8.3\n', netgen_path);
   fprintf(fid,'set TIX_LIBRARY=%s/lib/tcl8.2\n', netgen_path);
   fprintf(fid,'%s/ng431.exe %%*\n', netgen_path);
   fclose(fid);
end

disp('Netgen seems to have meshed your tank ok and written it to file!');
disp('..you just have to take your hats off to those guys at Johannes Kepler University, Linz,');
disp('what a good job.');

disp(['Now reading back data from file: ' meshfn])
[srf,vtx,fc,bc,simp,edg,mat_ind] = ng_read_mesh(meshfn);
disp([meshfn ' contains ' num2str(max(fc)) ' faces'])

%Translate to v3 notation!
tank_mdl.name = 'Tank model';
tank_mdl.nodes= vtx;
tank_mdl.elems= simp;
tank_mdl.boundary= srf;



%disp('Now I need some help finding which faces are electrodes')



% Plot wire frame equivalent of mesh model
%tetramesh(simp,vtx,'FaceColor','none','EdgeColor','cyan')%it doesn't work in matlab 5.3
if 0
set(gcf,'Name','Wire Mesh Model')
view(45,10)
hold on
trimesh(srf,vtx(:,1),vtx(:,2),vtx(:,3),'EdgeColor','blue')
title('Surface of body: blue, Volume of body: cyan')
hidden off
axis equal image; % Tightly fit square axes around plot
mshaxs = axis; % Save present axes for use with faces
pause(3)
end
% Select the electrodes
[elec,sels] = ng_tank_find_elec(srf,vtx,bc,centres);

%size(elec)

%[gnd_ind, electrodes, perm_sym, elec, protocol, no_pl] = get_model_elecs;


if size(elec,1) ~= nelec 
  error('That did''t work. Wrong number of electrodes!');
end

for i=1:nelec
    electrodes(i).z_contact= 1.0;  %well you can change that later!
    electrodes(i).nodes=     unique( elec(i,:) );
end

perm_sym='{n}';



tank_mdl.gnd_node=           1;
tank_mdl.electrode =         electrodes;
tank_mdl.misc.perm_sym =          perm_sym;

% save(fnstem,'tank_mdl');

