vi = fwd_solve(img);

img.elem_data(:) = 1;
vh = fwd_solve(img);

vi = add_noise(1e100, vi, vh);
