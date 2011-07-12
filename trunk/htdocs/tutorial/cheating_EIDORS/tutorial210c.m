% Test image for different noise
% $Id$

il_g= mk_common_model('c2c2',16);

num_tries=6;
for i= 1:num_tries 
   vi_n(:,i) = add_noise(2, vi, vh);
end

img = inv_solve( il_g, vh, vi_n );
img.calc_colours.greylev = .01;
img.show_slices.img_cols = num_tries;
show_slices( img );