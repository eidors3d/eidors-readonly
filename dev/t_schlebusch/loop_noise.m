%function loop_noise(r_real, cond_real, noiselevels)
run('C:\DATEN\eidors_3.7.1\startup.m')

%% Generate GREIT reconstruction
% GREIT Training
% Trainiere GREIT reconstruction matrix
opt.imgsz = [32 32];
%opt.distr = 3; %non-random, uniform
opt.distr = 2; %random, uniform
opt.Nsim = 1000;
%opt.target_size = 0.02; % Target size (frac of medium)
opt.target_size = [0.01 0.1]; % Target size (frac of medium)
opt.target_offset = [-0.1 0.1];
opt.noise_figure = 0.5; % Recommended NF=0.5;

fmdl_greit= ng_mk_cyl_models([0.3,0.15,0.01],[32,[0.1,0.2]],[0.01 0 0.002]); 
stim =  mk_stim_patterns(32,2,[0,1],[0,1], {'no_meas_current'}, 1);
fmdl_greit.stimulation = stim;
fmdl_greit = mdl_normalize(fmdl_greit, 0);
img_greit = mk_image(fmdl_greit,1);

imdl_greit = mk_GREIT_model(img_greit, 0.25, [], opt);

%% loop through algorithm comparison
%data = zeros(6,5);
%data(:,1) = [0.02:0.01:0.04];
%data(:,2) = 2;
volumina = [50:100:550]; % volume in ml
r_real = (volumina./((4/3*pi)*1000*1000)).^(1/3);
cond_real = [1:0.5:3];
noiselevels = [1e1 1e2 1e3 1e4];

r_real_rep = kron(r_real, ones(1,size(cond_real,2)));
cond_real_rep = repmat(cond_real', size(r_real,2), 1);
noiselevels_rep = kron(noiselevels, ones(1,size(r_real_rep,2)));
r_real_rep = repmat(r_real_rep', size(noiselevels',1), 1);
cond_real_rep = repmat(cond_real_rep, size(noiselevels',1), 1);


data = zeros(size(r_real_rep,1),6);
data(:,1) = r_real_rep;
data(:,2) = cond_real_rep;
data(:,3) = noiselevels_rep';

parallel.defaultClusterProfile('local');
c = parcluster();
job1 = createJob(c)

for r = 1:size(data,1)
    disp '------------------------------------------------------'
    disp(r);
    disp(data(r));
    disp '------------------------------------------------------'
    %run('C:\DATEN\eidors_3.7.1\startup.m')
    %path(path, 'C:\DATEN\eidors\dev\physics');
    %row = data(r,:);
    %[row(4), row(5), row(6)] = compare_noise(row(1), row(2), imdl_greit, row(3));
    %data(r,:) = row;
    %[data(r,1), data(r,2), data(r,3), data(r,4), data(r,5), data(r,6)] = compare_noise(data(r,1), data(r,2), imdl_greit, data(r,3));
    createTask(job1, @compare_noise, 6, {data(r,1), data(r,2), imdl_greit, data(r,3)})
end
job1.Tasks

%%
submit(job1)
wait(job1)
results = fetchOutputs(job1);
%disp(data)
save(strcat('loop_noise_',datestr(now,30),'.mat'))
%end