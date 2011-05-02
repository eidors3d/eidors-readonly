elec_per_plane = 8;
equal_spacing = 1;
elec_planes = [0.4,0.6];
elec_shape = [0.05]; % circular elecs,
maxh       = 0.08;   % mesh size
fmdl= mk_library_model({'adult_male','boundary'}, ...
       [elec_per_plane,equal_spacing, elec_planes], elec_shape, maxh);
subplot
subplot(211); show_fem(fmdl);
print_convert mk_library_model03a.png
view(2);
print_convert mk_library_model03b.png
