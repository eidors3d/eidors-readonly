%% Generate GREIT reconstruction
% GREIT Training
% Trainiere GREIT reconstruction matrix
opt.imgsz = [32 32];
opt.distr = 3; %non-random, uniform
opt.Nsim = 1000;
opt.target_size = 0.05; % Target size (frac of medium)
opt.noise_figure = 1.0; % Recommended NF=0.5;

fmdl_greit= ng_mk_cyl_models([0.3,0.15,0.005],[32,[0.1,0.2]],[0.01 0 0.002]); 
stim =  mk_stim_patterns(32,2,[0,1],[0,1], {'no_meas_current'}, 1);
fmdl_greit.stimulation = stim;
fmdl_greit = mdl_normalize(fmdl_greit, 0);
img_greit = mk_image(fmdl_greit,1);

imdl_greit = mk_GREIT_model(img_greit, 0.25, [], opt);

%% loop through algorithm comparison
%data = zeros(6,5);
%data(:,1) = [0.02:0.01:0.04];
%data(:,2) = 2;
r_real = [0.01:0.01:0.05];
cond_real = [0.5:0.5:3];

r_real_rep = kron(r_real, ones(1,size(cond_real,2)));
cond_real_rep = repmat(cond_real', size(r_real,2), 1);

data = zeros(size(r_real_rep,2),5);
data(:,1) = r_real_rep';
data(:,2) = cond_real_rep;

for r = 1:size(data,1)
    [data(r,3), data(r,4), data(r,5)] = compare(data(r,1), data(r,2), imdl_greit);
end
disp(data)