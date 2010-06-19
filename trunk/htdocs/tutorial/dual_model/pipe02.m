% Create pipe model $Id$

el= 'ellipsoid(0.3,0.2,1.5; 0,0,0.25; 0.1,0,0; 0,0.1,0)';
extra={'obj',['solid obj = ',el,';']};

[fmdl,mat_idx]= ng_mk_cyl_models(4,[n_elec,2],[0.2,0.5,0.04], extra);
fmdl.stimulation = stim;

img= mk_image(fmdl, 1);
vh = fwd_solve( img );
img.elem_data(mat_idx{2}) = 2;
vi = fwd_solve( img );

clf; subplot(121);
show_fem(img); view([0,0]);
print_convert('pipe02a.png','-density 100')

clf; plot( [vh.meas, 100*(vi.meas - vh.meas)] )
print_convert('pipe02b.png','-density 75')
