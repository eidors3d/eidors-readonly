


inv_crime=0;


%Stimulation - 16 elecs adjacent current/adjacent voltage
stim = mk_stim_patterns(16,1,'{ad}','{ad}');


%Make a uniform image with unit background
cond_bkg = 1;

fmdl= ng_mk_cyl_models([0,1,0.1],[16],[0.2,0,0.1]);
%fmdl= ng_mk_cyl_models(0,[16],[0.2,0,0.1]);

fmdl=fix_model(fmdl);
fmdl.stimulation = stim; %Add to model
fmdl.approx_type='tri3';

%Inverse model
imdl.solve='eidors_default';
%   imdl.hyperparameter=0.03;
imdl.RtR_prior='eidors_default';
imdl.jacobian_bkgnd.value=1;
imdl.reconst_type='difference';
imdl.fwd_model = fmdl;
imdl.name='built model';
imdl.type='inv_model';

%Notation of homog_img before
homog_img = mk_image(fmdl,cond_bkg);
homog_img.fwd_solve.get_all_meas=1;
homog_img.fwd_model.stimulation=stim;
% homog_img.fwd_model.M_tensor.a = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.b = ones(size(homog_img.elem_data,1),1);
% homog_img.fwd_model.M_tensor.rot = zeros(size(homog_img.elem_data,1),1);
fmdl = calc_closest_ellipse(fmdl);

data_hom = fwd_solve(homog_img);


%% Ptensor options
fmdl = calc_closest_ellipse(fmdl);
% Default options
% popts = [];
% popts.H0_type = 'ptensor';
% popts.use_hyper = 1;
% popts.neumann = 'disc';
% popts.ptensor_its = 100;
% popts.max_its = 100;
% popts.update_delta = 1;
% popts.inv_crime=inv_crime;
% popts.update_U0 = 1;
% popts.flexible = true;

%opts.neumann = 'disc';
popts.neumann = 'freespace';

popts.flexible = true;
popts.update_U0 = 1;


%% test scenarios
sepv=0.5:0.01:1;
radv=0.5;%logspace(log10(0.16), log10(0.25),2);
%%
% sepv=0.2
% radv =0.1646
% conv=1.5

cond_obj1 = 2;
cond_obj2 = 2;%1.5;

nsep = length(sepv);
nrad = length(radv);

img_pt = cell(nsep, nrad);
r_pt = img_pt;
er_f11 = img_pt;
img_eid_diff = img_pt;
img_eid_abs = img_pt;

t_pt = zeros(nsep, nrad);
t_diff = t_pt;
t_eabs = t_pt;
eres = t_pt;
cts_pt = t_pt;

ddiff = zeros(1,nsep);
mnorm = ddiff;
measindx = 140; % was biggest difference

for ii = 1:nsep
    
    fprintf('separation %f\n',sepv(ii));

    
    % Sim data
    sep = sepv(ii);
    rad = radv;
    sb1 = sprintf('solid ball1 = cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',sep,sep,sep,sep,rad);
    sb2 = sprintf('solid ball2 = cylinder(-%0.1f,-%0.1f,0;-%0.1f,-%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) and not cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) -maxh=0.1;',sep,sep,sep,sep,rad, sep,sep,sep,sep,rad);
    extra={'ball1','ball2',[sb1,sb2]};
    
    fmdl_t= ng_mk_cyl_models([0., 3, 0.05],[16],[0.2,0,0.05],extra);
    fmdl_t=fix_model(fmdl_t);
    fmdl_t.stimulation = stim; %Add to model
    fmdl_t.approx_type='tri3';
    sim_img{ii}= mk_image(fmdl_t, cond_bkg );
    sim_img{ii}.elem_data(fmdl_t.mat_idx{2})=cond_obj1;
    sim_img{ii}.elem_data(fmdl_t.mat_idx{3})=cond_obj2;
    sim_img{ii}.fwd_solve.get_all_meas = 1;
    %   figure; show_fem(img);
    pixel_group=[fmdl_t.mat_idx{2}];
    %pixel_group = [102,123,103,83,66,82]; %b2C
    %pixel_group = [327,364,328,292,259,291]; %c2C
    %pixel_group = 1:length(fmdl.elems(:,1));
    
    n_level = 50;
    
    data_sep{ii} = fwd_solve(sim_img{ii});
    
    
    sim_img{ii}.elem_data(fmdl_t.mat_idx{2})=cond_obj1;
    sim_img{ii}.elem_data(fmdl_t.mat_idx{3})=cond_bkg;
    data_1{ii} = fwd_solve(sim_img{ii});
    
    
    sim_img{ii}.elem_data(fmdl_t.mat_idx{2})=cond_bkg;
    sim_img{ii}.elem_data(fmdl_t.mat_idx{3})=cond_obj2;
    data_2{ii} = fwd_solve(sim_img{ii});
    
    vdiff{ii} = data_sep{ii}.volt(:,1) - data_1{ii}.volt(:,1) - data_2{ii}.volt(:,1);
    ddiff(ii) = norm(data_sep{ii}.meas(:) - (data_1{ii}.meas(:) + data_2{ii}.meas(:)));
    mnorm(ii) = norm(data_sep{ii}.meas);
            
end
%%
figure;
sim_img{1}.elem_data = vdiff{1};
show_fem(sim_img{1}, 1)
figure;
plot(sepv, ddiff./mnorm);

%% Contrast
conv = 0.1:0.01:10;
measindx = 69; % was max change
meas = zeros(1,length(conv));
radv = 0.2;

for ii=1:length(conv)
    fprintf('contrast %f\n',conv(ii));
    
    sep = 0.2;
    rad = radv;
    sb1 = sprintf('solid ball1 = cylinder(%0.1f,%0.1f,0;%0.1f,%0.1f,1;%0.1f) and orthobrick(-1,-1,0;1,1,0.05) -maxh=0.1;',sep,sep,sep,sep,rad);
    extra={'ball1',[sb1]};
    
    fmdl_t= ng_mk_cyl_models(0.,[16],[0.2,0,0.05],extra);
    fmdl_t=fix_model(fmdl_t);
    fmdl_t.stimulation = stim; %Add to model
    fmdl_t.approx_type='tri3';
    sim_img{ii}= mk_image(fmdl_t, cond_bkg );
    sim_img{ii}.elem_data(fmdl_t.mat_idx{2})=conv(ii);
    %   figure; show_fem(img);
    pixel_group=[fmdl_t.mat_idx{2}];
    %pixel_group = [102,123,103,83,66,82]; %b2C
    %pixel_group = [327,364,328,292,259,291]; %c2C
    %pixel_group = 1:length(fmdl.elems(:,1));
    
    n_level = 50;
    
    data_sat{ii} = fwd_solve(sim_img{ii});
    meas(ii) = data_sat{ii}.meas(measindx);
end


figure;
plot(conv, meas)
