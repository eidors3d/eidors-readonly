imgr = inv_solve(inv3d, vh, vi);

show_fem(imgr); view(-16,22); set(gca,'Projection','perspective')
print_convert surface_sim04a.png '-density 100'

% Set levels for z-intercepts of -0.4,-1.0,-1.6
levels= [inf,inf,-0.4,1,1; inf,inf,-1.0,2,1; inf,inf,-1.6,3,1];
show_slices(imgr,levels)
print_convert surface_sim04b.png '-density 150'
