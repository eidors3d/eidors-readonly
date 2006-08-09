function [tank_mdl,centres] = create_tank_mesh_ng( tank_radius, tank_height,CorR, log2_electrodes_per_plane, no_of_planes,first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem)
% USAGE: [tank_mdl,centres] = create_tank_mesh_ng( tank_radius, tank_height,CorR, log2_electrodes_per_plane, no_of_planes,first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem)
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
%
% $Id: create_tank_mesh_ng.m,v 1.12 2006-08-09 13:54:05 aadler Exp $

elecs_per_plane= 2^log2_electrodes_per_plane;


nelec = no_of_planes*elecs_per_plane;
if CorR=='C'
   electrode_radius=electrode_width; 
   if elecs_per_plane*electrode_radius*2 > 2*pi*tank_radius
      error('requested electrode width is larger than the plane')
   end
else
   if elecs_per_plane*electrode_width > 2*pi*tank_radius
      error('requested electrode width is larger than the plane')
   end
end


% Need some sanity checks here on the data!


geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
[fid,mess]=fopen(geofn,'w');

if CorR =='R'

   % Rectangular case 
   [fid,mess]=fopen(geofn,'w');

   write_header(fid,tank_height,tank_radius);

   kel = 0;
   for l=1:no_of_planes
      z = first_plane_starts + (l-1)*height_between_centres;
      for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/elecs_per_plane
          kel=kel+1;
          [x,y]=pol2cart(th,tank_radius);
          centres(kel,:)= [x,y,z]; % keep the centres
          dirn = [x,y,0];
          dirn = dirn ./ norm(dirn);
         
          % Write basic electrode shape
          writengcuboid(fid,sprintf('rod%d',kel),[x,y,z],dirn, ...
                        electrode_height,electrode_width,0.5*tank_radius);	 	 
          
          % Write surrounding box to force dense electrodes
%         writengcuboid(fid,sprintf('out_rod%d',kel),[x,y,z],dirn, ...
%                       1.05*electrode_height,1.1*electrode_width,0.5*tank_radius);	 	 
      end;
   end;
else

% Circular case
   write_header(fid,tank_height,tank_radius);


   kel = 0;
   for l=1:no_of_planes
      z = first_plane_starts + (l-1)*height_between_centres;
      for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/elecs_per_plane
         kel=kel+1;
         [x,y]=pol2cart(th+0.02*randn(1),tank_radius);  %jiggle to avoid netgen errors
         dirn = [x,y,0];
         centres(kel,:)= [x,y,z]; % keep the centres
         dirn = dirn ./ norm(dirn);
         writengcylrod(fid,sprintf('rod%d',kel),[x,y,z],dirn, ....
                       electrode_radius, 0.5*tank_radius);	 	 

%        writengcylrod(fid,sprintf('out_rod%d',kel),[x,y,z],dirn, ....
%                      1.1*electrode_radius, 0.5*tank_radius);	 	 
      end;
   end;

end % of circular case


   for k= 1:nelec
     fprintf(fid,'solid cyl%d = bigcyl    and rod%d; \n',k,k);
   end


   fprintf(fid,'tlo bigcyl;\n');
   for k=1:nelec
     fprintf(fid,'tlo cyl%d cyl -col=[1,0,0];\n ',k);
   end


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


function write_header(fid,tank_height,tank_radius);
   fprintf(fid,'#Automatically generated cylinder surface mesh with rectangular electrodes\n');
   fprintf(fid,'algebraic3d\n');
   fprintf(fid,'solid cyl=cylinder (0,0,0;0,0,%6.2f;%6.2f); \n', ...
           tank_height, tank_radius);
   fprintf(fid,['solid bigcyl= plane(0,0,0;0,0,-1)\n' ...
                'and  plane(0,0,%6.2f;0,0,1)\n' ...
                'and  cyl;\n'],tank_height);  


