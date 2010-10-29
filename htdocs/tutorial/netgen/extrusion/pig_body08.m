img = mk_image( fmdl, 1 );
img.elem_data( mat_idx{2} ) = 0.3;
vh= fwd_solve(img);

% Put a ball in the object center
targ= mk_c2f_circ_mapping(fmdl, [-2.8;-2.8;1;0.3]);
img.elem_data = img.elem_data + targ*.5;
subplot(234); show_fem(img); view(2);
title('Simulated change region');

show_fem(img); view(0,90);
print_convert pig_body08a.jpg


vi = fwd_solve(img);
vi = add_noise( 3, vi, vh );
plot([vh.meas, 100*(vi.meas - vh.meas)]);

print_convert pig_body08b.jpg
