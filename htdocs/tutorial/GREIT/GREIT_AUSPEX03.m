imdl = mk_common_model('c2c2',16);
imdl.fwd_model.stimulation =  ...
     mk_stim_patterns(16,1,[0,1],[0,1],{'rotate_meas','no_meas_current'},1);
img = mk_image(imdl);     vh = fwd_solve( img );
img.elem_data(290) = 1.1; vi = fwd_solve( img );

clf;subplot(221);
show_fem(img); title('Test Object');
print_convert GREIT_AUSPEX03a.png
