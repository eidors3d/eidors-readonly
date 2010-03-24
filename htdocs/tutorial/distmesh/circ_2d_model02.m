% Simple model with three electrodes $Id$
n_elecs= 14;
elec_width= 0.1;

hw= elec_width/2;
th = linspace(0,2*pi,n_elecs+1); th(end)=[];
for i=1:n_elecs;
   ti = th(i) + [hw;-hw];
   elec_pts{i} = [sin(ti),cos(ti)];
end

subplot(221);
fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.05] );
show_fem(fmdl);

subplot(222);
fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );
show_fem(fmdl);

print -dpng -r125 circ_2d_model02.png
