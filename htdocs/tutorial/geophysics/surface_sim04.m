imgr = inv_solve(inv3d, vh, vi);

show_fem(imgr); view(-16,22); set(gca,'Projection','perspective')
print_convert surface_sim04a.png '-density 100'

% Set levels for z-intercepts of -0.4,-1.0,-1.6
levels= [inf,inf,-0.4,1,1; inf,inf,-1.0,2,1; inf,inf,-1.6,3,1];
show_slices(imgr,levels)
print_convert surface_sim04b.png '-density 150'

% Try a different prior
inv3d.RtR_prior= @prior_noser;
inv3d.prior_use_fwd_not_rec = 1;
inv3d.hyperparameter.value = .3;

imgr = inv_solve(inv3d, vh, vi);

show_fem(imgr); view(-16,22); set(gca,'Projection','perspective')
print_convert surface_sim04c.png '-density 100'

show_slices(imgr,levels)
print_convert surface_sim04d.png '-density 150'
