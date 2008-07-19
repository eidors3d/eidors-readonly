% Create Model $Id: if_lung_injury01.m,v 1.1 2008-07-19 20:58:23 aadler Exp $

imdl= mk_common_model('c2c2',16);
% Reverse electrodes to give 'clinical' view (looking toward patient head)
imdl.fwd_model.electrode =  ...
   imdl.fwd_model.electrode([9:16,1:8]);
% Use normalized difference imaging
imdl.fwd_model.normalize_measurements=1;

subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.1,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,+1.1,'Dorsal');
set(hh,'HorizontalAlignment','Center');

print -r150 -dpng if_lung_injury01.png
