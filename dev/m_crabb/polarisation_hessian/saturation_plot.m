%Make an ivnerse model - standard foward model inside
imdl = mk_common_model('c2C',16);
fmdl = imdl.fwd_model; %Extract model

%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl.stimulation = stim; %Add to model

%Make an image and show image
img = mk_image(fmdl,1);
figure; show_fem(img,[0,1,3])

%Solve forward model
vh = fwd_solve(img);
figure; plot(vh.meas); 

%Pixk group of pixels to perturb
pixel_group = [327,364,328,292,259,291];

%Conductivity values
cond_vals = logspace(-2,2,100);

for jj=1:length(cond_vals)
    jj
    img.elem_data([pixel_group]) =cond_vals(jj);
    vi = fwd_solve(img); %Image with pixel perturbation
    norm_vi_diff(jj) = norm(vi.meas-vh.meas);
    
end




figure; semilogx(cond_vals,norm_vi_diff);
ylabel('|| Vi - V0||'); xlabel('Pertubation in ith pixel');