% generate noisy EIT voltage measurements 
% $ Id:  $

for ii = 1:length(imgs)
    imgs{ii}.elem_data = imgs{ii}.elem_data1;
    vh{ii} = fwd_solve(imgs{ii});
    imgs{ii}.elem_data = imgs{ii}.elem_data2;
    vi{ii} = fwd_solve(imgs{ii});    
end

% generate additive white Gaussian noise
NOISE_N_REALIZATIONS = 1E3;
NOISE_AMPLITUDE = 5E-3;
noise = NOISE_AMPLITUDE * randn(32*32, NOISE_N_REALIZATIONS);