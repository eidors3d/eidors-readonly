% Simulate obj $Id$

fmdl = ng_mk_cyl_models([2,1,0.08],[16,1],[0.05]); 
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
imgs= mk_image( fmdl, 1);

show_fem(imgs);
print_convert('GREIT_test_params01a.png','-density 60');

view(0,0)
xlim([-.4,.4])
zlim(1+[-.4,.4])
print_convert('GREIT_test_params01b.png','-density 60');
