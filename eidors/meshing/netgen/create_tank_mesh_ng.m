function [tank_mdl,centres] = create_tank_mesh_ng( ...
            tank_radius, tank_height,CorR, log2_electrodes_per_plane, ...
            no_of_planes,first_plane_starts, height_between_centres,  ...
            electrode_width, electrode_height,fnstem, elec_mesh_density )
%USAGE: [tank_mdl,centres] = create_tank_mesh_ng( ...
%           tank_radius, tank_height,CorR, log2_electrodes_per_plane, ...
%           no_of_planes,first_plane_starts, height_between_centres,  ...
%           electrode_width, electrode_height,fnstem, elec_mesh_density )
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
% elec_mesh_density  - force enhanced mesh density near electrodes
%           default is to not enhance mesh density
%
%
%
% Function to generate tank model
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
% Revised for new version 3.0 structure WRBL 05/12/2005
%Made in to function WRBL 6/5/2005
%
% $Id$

elecs_per_plane= 2^log2_electrodes_per_plane;

% meshsize around electrodes
if nargin<11
   elec_mesh_density= 0; % points on outsize of electrode
end

nelec = no_of_planes*elecs_per_plane;

% Need some sanity checks here on the data!

geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
sizefn= [fnstem,'.msz'];
fid=fopen(geofn,'w');

if elec_mesh_density > 0
   fsf=fopen(sizefn,'w');
end


write_header(fid,tank_height,tank_radius);

switch CorR
   case 'R'
      if elecs_per_plane*electrode_width > 2*pi*tank_radius
         error('requested electrode width is larger than the plane')
      end
      [centres] = write_rectangular_elecs(fid, ...
               elecs_per_plane, tank_radius, tank_height, ...
               no_of_planes,first_plane_starts, height_between_centres,  ...
               electrode_width, electrode_height,fnstem, elec_mesh_density );
   case 'C'
      electrode_radius=electrode_width; 
      if elecs_per_plane*electrode_radius*2 > 2*pi*tank_radius
         error('requested electrode width is larger than the plane')
      end
      [centres] = write_circular_elecs(fid, ...
               elecs_per_plane, tank_radius, tank_height, ...
               no_of_planes,first_plane_starts, height_between_centres,  ...
               electrode_radius, fnstem, elec_mesh_density );
   otherwise;
      error('Specified electrode must be C or R');
end

   for k= 1:nelec
     fprintf(fid,'solid cyl%d = bigcyl    and rod%d; \n',k,k);
   end


   fprintf(fid,'tlo bigcyl;\n');
   for k=1:nelec
     fprintf(fid,'tlo cyl%d cyl -col=[1,0,0];\n ',k);
   end


fclose(fid);
if elec_mesh_density>0
   fclose(fsf);
end

% Now call Netgen in batchmode to mesh this CSG file
disp('Calling Netgen. Please wait.....');
if elec_mesh_density>0
   call_netgen( geofn, meshfn, sizefn);
else
   call_netgen( geofn, meshfn);
end


disp('Netgen seems to have meshed your tank ok and written it to file!');
disp('..you just have to take your hats off to those guys at Johannes Kepler University, Linz, what a good job.');
disp('Please note that netgen has some strange limitations that we have not been able to fully debug. To get around this, we have added a small random jitter to the config files, which randomly solves the problem. Thus, you may need to try to re-run your program a few times, if you are getting strange errors');

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


function [centres] = write_rectangular_elecs(fid, ...
            elecs_per_plane, tank_radius, tank_height, ...
            no_of_planes,first_plane_starts, height_between_centres,  ...
            electrode_width, electrode_height,fnstem, elec_mesh_density )

   if elec_mesh_density>0
      mesh_density= (electrode_height+electrode_width)/elec_mesh_density;
      %Density has some very funny limits. It sometimes just
      % breaks. Try adding random jitter so that it will
      % at least work sometimes - or after multiple attempts
      mesh_density= mesh_density*(1+randn(1)*.01);
   
      % add 2 points because we duplicate the points at each corner
      vert_points= ceil( elec_mesh_density*electrode_height / ...
                         2 / (electrode_height + electrode_width) ) + 2;
      vert_sides= linspace(-electrode_height/2, ...
                            electrode_height/2, vert_points );

      horz_points= ceil( elec_mesh_density*electrode_width / ...
                         2 / (electrode_height + electrode_width) ) + 2;
      horz_sides= linspace(-electrode_width/2, ...
                            electrode_width/2, horz_points );

      fprintf(fsf,'%d\n', no_of_planes*elecs_per_plane* ...
               (vert_points + horz_points)*2 );

      yy= [vert_sides,vert_sides, ...
           electrode_height/2*ones(1,horz_points), ...
          -electrode_height/2*ones(1,horz_points)]';
      zz= [electrode_width/2*ones(1,vert_points), ...
          -electrode_width/2*ones(1,vert_points), ...
           horz_sides,horz_sides]';
      th_col= 0*yy+1;
      xyz= [0*th_col,yy,zz];
   end


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
          if elec_mesh_density>0
             %write to meshsize file
             xyz1= th_col*[x,y,z] + xyz* ...
                       [cos(th),sin(th),0; ...
                       -sin(th),cos(th),0; ...
                        0      ,0      ,1];
             fprintf(fsf,'%f %f %f %.3f\n',[xyz1,th_col*mesh_density]');
%            plot3(xyz1(:,1),xyz1(:,2),xyz1(:,3),'.') % View to debug
          end

      end;
   end;

function [centres] = write_circular_elecs(fid, ...
            elecs_per_plane, tank_radius, tank_height, ...
            no_of_planes,first_plane_starts, height_between_centres,  ...
            electrode_radius, fnstem, elec_mesh_density )

   if elec_mesh_density>0
      fprintf(fsf,'%d\n', no_of_planes*elecs_per_plane*elec_mesh_density);
      mesh_density= pi*electrode_radius/elec_mesh_density;
      %Density has some very funny limits. It sometimes just
      % breaks. Try adding random jitter so that it will
      % at least work sometimes
      mesh_density= mesh_density*(1+randn(1)*.01);
   end

   
   theta= linspace(0,2*pi, elec_mesh_density+1)'; theta(1)=[];
   th_col= ones(size(theta));
   [xyz]= [0*th_col, ...
           electrode_radius*sin(theta), ...
           electrode_radius*cos(theta)];

   xyz2= [];
   jiggle= 0.00; %jiggle to avoid netgen errors
   kel = 0;
   for l=1:no_of_planes
      % test effect of this

      z = first_plane_starts + (l-1)*height_between_centres;
      for th =0:2*pi/elecs_per_plane:2*pi*(elecs_per_plane-1)/elecs_per_plane
         kel=kel+1;
%        modify a single electrode to see the effect
%        if kel ==elecs_per_plane; electrode_radius=electrode_radius*0.7; end
         [x,y]=pol2cart(th+jiggle,tank_radius); 
         dirn = [x,y,0];
         centres(kel,:)= [x,y,z]; % keep the centres
         dirn = dirn ./ norm(dirn);
         writengcylrod(fid,sprintf('rod%d',kel),[x,y,z],dirn, ....
                       electrode_radius, 0.2*tank_radius);

         if elec_mesh_density>0
            %write to meshsize file
            xyz1= th_col*[x,y,z] + xyz* ...
                      [cos(th),sin(th),0; ...
                      -sin(th),cos(th),0; ...
                       0      ,0      ,1];
            fprintf(fsf,'%f %f %f %.3f\n',[xyz1,th_col*mesh_density]');
         end

      end;
   end;

function writengcuboid(fid,name,c, dirn,h,w,d)
% writes the specification for a netgen cuboid on fid, named name, centerd on c,
% in the direction given by vector dirn, height  h width w and depth d
% direction is in the xy plane
dirnp = [-dirn(2),dirn(1),0];
bl = c - (d/2)* dirn + (w/2)* dirnp -[0,0,h/2];
tr =c + (d/2)* dirn - (w/2)* dirnp +[0,0,h/2];
fprintf(fid,'solid %s  =plane (%6.3f,%6.3f,%6.3f;0, 0, -1 )\n ', ...
        name,bl(1),bl(2),bl(3));
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
        bl(1),bl(2),bl(3),-dirn(1),-dirn(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
        bl(1),bl(2),bl(3),dirnp(1),dirnp(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;0, 0, 1  )\n ', ...
        tr(1),tr(2),tr(3));
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
        tr(1),tr(2),tr(3),dirn(1),dirn(2),0);
fprintf(fid,'        and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  );\n ', ...
        tr(1),tr(2),tr(3),-dirnp(1),-dirnp(2),0);

function writengcylrod(fid,name,c, d,rd,ln)
% writes the specification for a netgen cylindrical rod on fid, named name, centerd on c,
% in the direction given by vector d, radius rd  lenght ln
% direction is in the xy plane
 % the direction vector
dirn = d.*(ln/(2*norm(d)));   % normalize
inpt = c - dirn.*(ln/2);
outpt =c + dirn.*(ln/2);


fprintf(fid,'solid %s  =plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
      name , inpt(1),inpt(2),inpt(3),-dirn(1),-dirn(2),-dirn(3));
fprintf(fid,'       and plane(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f  )\n ', ...
      outpt(1),outpt(2),outpt(3),dirn(1),dirn(2),dirn(3));
fprintf(fid,'       and cylinder(%6.3f,%6.3f,%6.3f;%6.3f,%6.3f,%6.3f; %6.3f  );\n ', ...
      inpt(1),inpt(2),inpt(3),outpt(1),outpt(2),outpt(3), rd);
