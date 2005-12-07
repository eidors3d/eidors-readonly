% Script to generate tank model interactively
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
% Revised for new version 3.0 structure BL 05/12/2005
% $Id: create_tank_mesh_ng_interactive.m,v 1.3 2005-12-07 22:04:04 aadler Exp $
disp('EIDORS 3D tank generation script');
disp('You have to answer lots of tedious questions');
disp('Someone please make a GUI!');

while 1 % loop until user is happy with selection
CorR = upper(input('Circular electrodes or rectangular? [C/R]','s'));
tank_radius = input('Tank radius? ');
tank_height = input('Tank height? ');
disp('Number of electrodes on each plane is 2^k? ');
elec_per_plane_index=input('Input the index k? ');
elecs_per_plane= 2^elec_per_plane_index;
log2_electrodes_per_plane=elec_per_plane_index;
disp(sprintf('Thats %d electrodes per plane? ',elecs_per_plane));
no_of_planes = input('Number of planes? ');


nelec = no_of_planes*elecs_per_plane;
disp(sprintf('Ok so you have %d electrodes in total? ',nelec));
if nelec < 9
  disp('That''s not very many!');
elseif nelec <17
  disp('Come on this is 3D EIT. It''s not the 1980s you know!');
elseif nelec <65
  disp('Ok sounds like enough!');
else
  disp('Awesome!');
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

% Print user's selections. Ask if they're happy
disp('----------------------------------------');
disp('SUMMARY OF YOUR SELECTIONS:');
disp('----------------------------------------');
fprintf('Tank radius?                        = %f\n',tank_radius);
fprintf('Tank height?                        = %f\n',tank_height);
fprintf('Number of electrodes on each plane  = %d\n',elecs_per_plane);
fprintf('Number of planes?                   = %d\n',no_of_planes);
fprintf('Height of centre of first plane?    = %d\n',first_plane_starts); 
fprintf('Height between centres?             = %d\n',height_between_centres); 
if CorR=='C'
disp('Circular Electrodes');
fprintf('Electrode radius?                   = %f\n',electrode_radius);
else
disp('Rectangular Electrodes');
fprintf('Electrode height?                   = %f\n',electrode_height);
fprintf('Electrode width?                    = %f\n',electrode_width );
end
fprintf('Base filename for .geo and .vol file= %s\n', fnstem);
disp('----------------------------------------');
OK= input('Are these inputs correct? [Y/n]','s');
if upper(OK)~='N' 
    break;
end

end % repeat back to while

geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
[fid,mess]=fopen(geofn,'w');

[tank_mdl,centres] = create_tank_mesh_ng( tank_radius, tank_height, CorR,...
log2_electrodes_per_plane, no_of_planes,...
first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem);
save(fnstem,'tank_mdl')

