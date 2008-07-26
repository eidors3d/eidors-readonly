% Dual Partial $Id$

imdl= mk_common_model('d2c2',64);
fmdl= imdl.fwd_model;
% Only keep 11 electrodes on top
fmdl.electrode = fmdl.electrode([60:64,1:6]);
% New stimulation pattern with 11 electrodes
fmdl.stimulation = mk_stim_patterns(11,1,'{ad}','{ad}',{},1);
% Remove meas_select - it was created for 64 electrodes
fmdl = rmfield(fmdl,'meas_select'); 
subplot(121)
show_fem(fmdl); axis square

% Crop model
[cmdl,c2f_idx]= crop_model(fmdl, inline('(y-0.25)-abs(x)>0','x','y','z'));
subplot(122)
show_fem(cmdl);
axis(1.05*[-1 1 -1 1]); axis square;

print -r125 -dpng dual_partial2d01a.png;
