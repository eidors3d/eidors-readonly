% load some lung data
load iirc_data_2006;
vi = v_rotate(:,50);
vh = v_reference;

img1= inv_solve(imdl_e,vh,vi);
subplot(121)
show_slices(img1);

img2= inv_solve(imdl_n,vh,vi);
subplot(122)
show_slices(img2);


print -r125 -dpng nodal_jacobian_solver02.png;

