% Call Netgen $Id: meshsize_param01.m,v 1.1 2008-06-11 18:42:39 aadler Exp $ 
electrodes_per_plane= 16;
number_of_planes= 1;
movement_pattern='radial_turn';
tank_height= 10;
finelevel=''; %finelevel='-fine'; Other option
refine_electrodes=0;
fno_max=1;
maxh= ''; %don't use max h
move_the_ball % CALL NETGEN
