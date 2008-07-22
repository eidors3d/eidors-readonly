% Create Model $Id: if_lung_injury01.m,v 1.4 2008-07-22 15:52:41 aadler Exp $

imdl= mk_common_model('c2c2',16);
% Reverse electrodes to give 'clinical' view (looking toward patient head)
imdl.fwd_model.electrode =  ...
   imdl.fwd_model.electrode([9:-1:1,16:-1:10]);
% Use normalized difference imaging
imdl.fwd_model.normalize_measurements=1;

subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.15,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,+1.15,'Ventral');
set(hh,'HorizontalAlignment','Center');

print -r150 -dpng if_lung_injury01.png
