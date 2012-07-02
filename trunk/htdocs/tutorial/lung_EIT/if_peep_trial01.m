% Create Model $Id$

fmdl = mk_library_model('pig_23kg_16el');
[fmdl.stimulation fmdl.meas_select] = mk_stim_patterns(16,1,'{ad}','{ad}');
fmdl = mdl_normalize(fmdl, 1);  % Use normalized difference imaging
opt.noise_figure = 0.5; opt.imgsz = [64 64];
imdl = mk_GREIT_model(fmdl, 0.25, [], opt);

% subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.1,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,-1.15,'Ventral');
set(hh,'HorizontalAlignment','Center');

print_convert if_peep_trial01.png
