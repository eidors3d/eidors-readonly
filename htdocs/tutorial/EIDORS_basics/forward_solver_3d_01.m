nelec= 7;
fmdl= ng_mk_cyl_models([1,0.3,0.05],[nelec,0.2],[0.1,0.05,0.02]);

stim = mk_stim_patterns(nelec,1,[0,1],[0,1],{'meas_current'},1);
fmdl.stimulation = stim;
conduct = 1;
img = mk_image( fmdl, conduct ); 
img.fwd_solve.get_all_meas = 1;
vh = fwd_solve(img);

img_v = rmfield(img, 'elem_data');
img_v.node_data = vh.volt(:,1);

subplot(231); show_fem(img);
subplot(232): plot(vh.meas)
subplot(233);
img_v.calc_colours.npoints=128;
show_slices(img_v,[inf,inf,0.2])


extra={'ball','solid ball = sphere(0,0,0.5;0.1);'};
[fmdl,mat_idx]= ng_mk_cyl_models([1,0.3,0.05],[nelec,0.2],[0.1,0.05,0.02],extra);
fmdl.stimulation = stim;

img= mk_image(fmdl, conduct);
img.fwd_solve.get_all_meas = 1;
vh2 = fwd_solve(img);

img.elem_data(mat_idx{2}) = 2;
vi2 = fwd_solve(img);



subplot(234); show_fem(img);
subplot(235): plot([vi2.meas, 100*(vh2.meas-vi2.meas)])


img_v = rmfield(img, 'elem_data');
img_v.node_data = vi2.volt(:,1) - vh2.volt(:,1);

subplot(236);
img_v.calc_colours.npoints=128;
show_slices(img_v,[inf,inf,0.2])
