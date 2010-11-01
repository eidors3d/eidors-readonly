load CT1

[fmdl, mat_idx] = ng_mk_extruded_model({2,{thorax,rlung,llung},[4,50],.1},[16,1,1],[.1,0]);
[stim,meas_sel] = mk_stim_patterns(16,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdl.stimulation = stim;

fmdl.nodes = fmdl.nodes*diag([ 1, 1,1]); % Flip x,y axis to match medical direction
img=mk_image(fmdl,1);
 img.elem_data(mat_idx{2})= 0.3;
 img.elem_data(mat_idx{3})= 0.3;

clf; show_fem(img); view(0,70);
