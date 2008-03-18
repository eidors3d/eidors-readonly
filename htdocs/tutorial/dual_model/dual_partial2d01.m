% Dual Partial $Id: dual_partial2d01.m,v 1.1 2008-03-18 16:49:16 aadler Exp $

imdl= mk_common_model('d2c2',64);
fmdl= imdl.fwd_model;
% Only keep 7 electrodes on top
fmdl.electrode = fmdl.electrode([62:64,1:4]);
% New stimulation pattern with 7 electrodes
fmdl.stimulation = mk_stim_patterns(7,1,'{ad}','{ad}',{},1);
subplot(121)
show_fem(fmdl); axis square

% Crop model
cmdl= crop_model(fmdl, inline('(y-0.45)-abs(x)>0','x','y','z'));
subplot(122)
show_fem(cmdl);
axis(1.05*[-1 1 -1 1]); axis square;

print -r125 -dpng dual_partial2d01a.png;
