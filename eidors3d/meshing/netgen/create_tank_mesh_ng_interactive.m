% Script to generate tank model interactively
% Bill Lionheart 23/01/2005 (somewhere over Siberia)
% Part of EIDORS 3D
% Revised for new version 3.0 structure BL 05/12/2005
disp('EIDORS 3D tank generation script');
disp('You have to answer lots of tedious questions');
disp('Someone please make a GUI!');


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

geofn= [fnstem,'.geo'];
meshfn= [fnstem,'.vol'];
[fid,mess]=fopen(geofn,'w');

[tank_mdl,centres] = create_tank_mesh_ng( tank_radius, tank_height, CorR,... log2_electrodes_per_plane, no_of_planes,...
first_plane_starts, height_between_centres, electrode_width, electrode_height,fnstem)
save(fnstem,'tank_mdl')

