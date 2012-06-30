vi = fwd_solve(img);

img.elem_data(:) = 1;
vh = fwd_solve(img);

vi = add_noise(5, vi, vh); %SNR=5
