function tank_mdl = create_tank_mesh_ng_centres( tank_radius, tank_height,CorR, centres, electrode_width, electrode_height,fnstem)
% USAGE: tank_mdl = create_tank_mesh_ng_centres( tank_radius, tank_height,CorR, centres,electrode_width, electrode_height, fnstem) 
%
% Parameters  tank_radius, tank_height, 
%   CorR= 'C' for circular 'R' for rectangular
% centres coordinates of centres of electrodes
% electrode_width, electrode_height , the width is just radius if 'R'
%  fnstem the file name used for saving netgen files
% Function to generate tank model
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
% Revised for new version 3.0 structure WRBL 05/12/2005
%Made in to function WRBL 6/5/2005
%
% $Id: create_tank_mesh_ng_centres.m,v 1.1 2006-08-03 06:56:43 billlion Exp $




if CorR=='C'
   electrode_radius=electrode_width; 
end


% Need some sanity checks here on the data!


geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
[fid,mess]=fopen(geofn,'w');

nelec = size(centres,1);

if CorR =='R'

% Rectangular case 
[fid,mess]=fopen(geofn,'w');
fprintf(fid,'#Automatically generated cylinder surface mesh with rectangular electrodes\n');
fprintf(fid,'algebraic3d\n');
fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,%6.2f;%6.2f); \n',tank_height, ...
	tank_radius);
fprintf(fid, 'solid bigcyl= plane(0,0,0;0,0,-1)\n  and  plane(0,0,%6.2f;0,0,1)\n   and  cyl;\n',tank_height);  
for kel = 1:nelec 
    
    x=centres(kel,1);
    y=centres(kel,2);
    z=centres(kel,3);
   
    dirn = [x,y,0];
    dirn = dirn ./ norm(dirn);
    writengcuboid(fid,sprintf('rod%d',kel),[x,y,z],dirn,electrode_height,electrode_width,0.5*tank_radius);	 	 
%    fprintf(fid,'solid cyl%d = bigcyl and  cylinder %(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f);\n',kcyl,-x,-y,z,x,y,z,electrode_radius );
end



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


for kel=1:nelec
    x=centres(kel,1);
    y=centres(kel,2);
    z=centres(kel,3);
    dirn = [x,y,0];
    dirn = dirn ./ norm(dirn);
    writengcylrod(fid,sprintf('rod%d',kel),[x,y,z],dirn,electrode_radius, 0.3*tank_radius);	 	 
%    fprintf(fid,'solid cyl%d = bigcyl and  cylinder %(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f;%6.3f);\n',kcyl,-x,-y,z,x,y,z,electrode_radius );
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
call_netgen( geofn, meshfn );

disp('Netgen seems to have meshed your tank ok and written it to file!');
disp('..you just have to take your hats off to those guys at Johannes Kepler University, Linz,');
disp('what a good job.');

disp(['Now reading back data from file: ' meshfn])
tank_mdl= ng_mk_fwd_model( meshfn, centres, ...
             'Netgen based cylindrical tank model', [] );


