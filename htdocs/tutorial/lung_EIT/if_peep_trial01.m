% Create Model $Id$

imdl= mk_common_model('c2c2',16);
% Use normalized difference imaging
imdl.fwd_model.normalize_measurements=1;

subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.1,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,+1.1,'Ventral');
set(hh,'HorizontalAlignment','Center');

print -r150 -dpng if_peep_trial01.png
