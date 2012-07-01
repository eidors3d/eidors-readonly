% Create Model $Id$

imdl= mk_common_model('c2c2',16);
% Reverse electrodes to give 'clinical' view (looking toward patient head)
imdl.fwd_model.electrode =  ...
   imdl.fwd_model.electrode([9:-1:1,16:-1:10]);
% Use normalized difference imaging
imdl.fwd_model = mdl_normalize(imdl.fwd_model, 1);

subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.15,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,+1.15,'Ventral');
set(hh,'HorizontalAlignment','Center');

print_convert if_lung_injury01.png 
