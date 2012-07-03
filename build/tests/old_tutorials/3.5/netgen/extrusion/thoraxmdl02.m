f = 100;
[fmdl, mat_idx] = ng_mk_extruded_model({2,{thorax/f,rlung/f,llung/f},[4,50],.1},[16,1.00,1],[.1,0,.05]);
fmdl.nodes = fmdl.nodes*f;
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;

img=mk_image(fmdl,1);
img.elem_data(mat_idx{2})= 0.3; % rlung
img.elem_data(mat_idx{3})= 0.3; % llung

clf; show_fem(img); view(0,70);
print_convert thoraxmdl02a.jpg
