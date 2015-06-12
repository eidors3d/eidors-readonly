clf
n_elecs = 16;

% Elliptic model
[fmdle,midx] = ng_mk_ellip_models([1, 1.14,1,0.15] ,[n_elecs,0.5],[0.05]);
[stim,msel] =  mk_stim_patterns(n_elecs,1,[0,1],[0,1],{'no_meas_current'}, 1);
fmdle.stimulation =  stim;       
fmdle.meas_select =  msel;       
fmdle = mdl_normalize(fmdle, 1);

% GREIT Ellip - circ objects
opt.distr = 0; % central
opt.noise_figure = 0.5;
imdl = mk_GREIT_model(mk_image(fmdle,1), 0.25, [], opt);

vh = mean(vv,2);        % reference is average
vi = vv(:,[45,70,173]); %3 inspirations

img = inv_solve(imdl,vh,vi);

img.show_slices.img_cols = 3;
img.show_slices.sep      = 2;
img.calc_colours.ref_level=0;
show_slices(img);

print_convert neonate_ex04a.png
