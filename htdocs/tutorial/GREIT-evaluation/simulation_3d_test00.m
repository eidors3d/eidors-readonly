% Mk 3D Netgen Model
electrodes_per_plane=16;
number_of_planes= 1;
tank_height= 30;
finelevel=' '; refine_electrodes= 0;
fno_max=0;
maxh= '-maxh=0.75'; % 268314 elems
maxh= '-maxh=1.5'; % 43164 elems
move_the_ball;
save ng_tank fmdl
