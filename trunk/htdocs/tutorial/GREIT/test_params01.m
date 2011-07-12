% create 3D forward model
fmdl = ng_mk_cyl_models([2,1,0.08],[8,0.8,1.2],[0.05]); 
fmdl.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{},1);
imgs= mk_image( fmdl, 1);

show_fem(imgs);
print_convert test_params01a.png '-density 50'
