fmdl.stimulation = stim_meas_list([1,num_elecs(fmdl),1,2]);

img= mk_image(fmdl,1);
img.fwd_solve.get_all_meas = true;
vh=fwd_solve(img);

imgh = rmfield(img,'elem_data');
imgh.node_data = vh.volt;
show_slices(imgh,1); %center slice

print_convert internal_electrodes03a.jpg
