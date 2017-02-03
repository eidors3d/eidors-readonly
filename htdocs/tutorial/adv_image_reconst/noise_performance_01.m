% create forward models of the human thorax 
% $Id$

% generate base model for forward solving
fmdl = mk_library_model('adult_male_32el_lungs');
imgBasic = mk_image(fmdl, 0.2);     % back ground conductivity
imgBasic.elem_data(fmdl.mat_idx{2}) = 0.13;   % left lung
imgBasic.elem_data(fmdl.mat_idx{3}) = 0.13;   % right lung
imgBasic.elem_data1 = imgBasic.elem_data;
% now generate a conductivity change
imgBasic.elem_data2 = imgBasic.elem_data;
imgBasic.elem_data2(fmdl.mat_idx{2}) = 0.1 * imgBasic.elem_data2(fmdl.mat_idx{2});
imgBasic.elem_data2(fmdl.mat_idx{3}) = 0.05 * imgBasic.elem_data2(fmdl.mat_idx{3}); 

% generate base model for inverse model
rmdlBasic = mk_library_model('adult_male_32el');

% 16 electrodes, equidistantly spaced, skip 0 (adjacent) stim/meas pattern
imgs{1} = imgBasic;
rmdls{1} = rmdlBasic;
imgs{1}.fwd_model.electrode(2:2:end) = []; 
imgs{1}.fwd_model.name = '16 elecs, skip 0';
rmdls{1}.electrode(2:2:end) = [];
imgs{1}.fwd_model.stimulation = mk_stim_patterns(16,1,[0,1],[0,1],{'no_rotate_meas','no_meas_current'});
rmdls{1}.stimulation = imgs{1}.fwd_model.stimulation;

% 16 electrodes, equidistantly spaced, stim/meas pattern with skip=5
imgs{2} = imgBasic;
rmdls{2} = rmdlBasic;
imgs{2}.fwd_model.electrode(2:2:end) = []; 
imgs{2}.fwd_model.name = '16 elecs, skip 5';
rmdls{2}.electrode(2:2:end) = [];
imgs{2}.fwd_model.stimulation = mk_stim_patterns(16,1,[0,1+5],[0,1+5],{'no_rotate_meas','no_meas_current'});
rmdls{2}.stimulation = imgs{2}.fwd_model.stimulation;

% 32 electrodes, equidistantly spaced, skip 0 (adjacent) stim/meas pattern
imgs{3} = imgBasic;
rmdls{3} = rmdlBasic;
imgs{3}.fwd_model.name = '32 elecs, skip 0';
imgs{3}.fwd_model.stimulation = mk_stim_patterns(32,1,[0,1],[0,1],{'no_rotate_meas','no_meas_current'});
rmdls{3}.stimulation = imgs{3}.fwd_model.stimulation;

% 24 electrodes, more densly spaced ventrally, stim/meas pattern with skip=9
imgs{4} = imgBasic;
rmdls{4} = rmdlBasic;
imgs{4}.fwd_model.electrode(9:2:24) = [];  % remove every 2nd elec on the back
rmdls{4}.electrode(9:2:24) = [];
imgs{4}.fwd_model.name = '24 elecs, skip 9';
imgs{4}.fwd_model.stimulation = mk_stim_patterns(24,1,[0,1+9],[0,1+9],{'no_rotate_meas','no_meas_current'});
rmdls{4}.stimulation = imgs{4}.fwd_model.stimulation;

% now plot every model 
clf;
for ii = 1:length(imgs)
    subplot(1,4,ii);
    h = show_fem(imgs{ii});
    set(h, 'edgecolor', 'none');
    title(imgs{ii}.fwd_model.name);
    view(2);    
    set(gca, 'ytick', [], 'xtick', []);
    xlabel('Dorsal');
    ylabel('Right');
end

print_convert np_models.png '-density 100'