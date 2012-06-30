% Call Netgen $Id$ 
stim = mk_stim_patterns(16,1,'{ad}','{ad}',{},1);
fmdl = ng_mk_cyl_models([10,15],[16,5],[0.5]);
fmdl.stimulation= stim;
img = mk_image(fmdl,1);
vh = fwd_solve(img);

extra={'ball','solid ball = sphere(7.5,0,5;2);'};
fmdl = ng_mk_cyl_models([10,15],[16,5],[0.5],extra);
fmdl.stimulation= stim;
img = mk_image(fmdl,1);
img.elem_data(fmdl.mat_idx{2}) = 1.1;
vi = fwd_solve(img);
