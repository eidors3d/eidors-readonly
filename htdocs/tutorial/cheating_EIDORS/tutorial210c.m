% Test image for different noise
% $Id: tutorial210c.m,v 1.1 2007-06-15 18:12:03 aadler Exp $

il_g= mk_common_model('c2c');

num_tries=6;
vi_n(1:num_tries) = vi;
for i= 1:num_tries % matlab can't vectorize
   noise = .002*randn(size(vi.meas));
   vi_n(i).meas = vi_n(i).meas + noise;
end

levels= [0,0,0,1,1];
show_slices( inv_solve( il_g, vi_n, vh ), levels);
