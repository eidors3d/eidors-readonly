% Create 2D model of a tunnel $Id$ 

extra={'ballorcyl',
  [ 'solid ball = sphere(0,0,0;1) -maxh=0.25;', ...
    'solid cylc = cylinder(0,0,0;0,0,1;6) -maxh=0.5;', ...
    'solid cylo = cylc and orthobrick(-10,-10,0;10,10,1);', ...
    'solid ballorcyl = ball or cylo;']};
extra={'ball', 'solid ball = sphere(0,0,0;1) -maxh=0.5;'};
cmdl= ng_mk_cyl_models([0,15,3],[0],[0.1,0,0.05],extra);

show_fem(cmdl);                      print_convert tunnelsim05a.png
show_fem(cmdl); axis(2*[-1,1,-1,1]); print_convert tunnelsim05b.png
