% Create Model $Id: if_peep_trial01.m,v 1.1 2008-07-19 16:42:22 aadler Exp $

imdl= mk_common_model('c2c2',16);
% Reverse electrodes to give 'clinical' view (looking toward patient head)
imdl.fwd_model.electrode =  ...
   imdl.fwd_model.electrode([1,16:-1:2]); 

subplot(211);
show_fem(imdl.fwd_model, [0,1,0]);
axis equal; axis off

hh=text(-1.1,0,'Right');
set(hh,'Rotation',90,'HorizontalAlignment','Center');
hh=text(0,+1.1,'Top');
set(hh,'HorizontalAlignment','Center');

print -r100 -dpng if_peep_trial01.png
