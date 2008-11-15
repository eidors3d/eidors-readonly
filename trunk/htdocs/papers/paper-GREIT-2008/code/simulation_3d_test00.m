% Mk 3D Netgen Model $Id$
electrodes_per_plane=16; number_of_planes= 1;
tank_height= 30;
refine_electrodes= 0; finelevel=' ';
 maxh= '-maxh=2.0'; %  29k elems
 maxh= '-maxh=1.5'; %  43k elems
%maxh= '-maxh=1.0'; % 124k elems
fno_max=0;
move_the_ball;
save ng_tank fmdl
