% Simulate object moving $Id: dual_partial2d02.m,v 1.1 2008-03-18 16:49:17 aadler Exp $

% simulation model
imdl= mk_common_model('d2c2',64);
smdl= imdl.fwd_model;
% Only keep 7 electrodes on top 
smdl.electrode = smdl.electrode([62:64,1:4]);
smdl.stimulation = mk_stim_patterns(7,1,'{ad}','{ad}',{},1);

[vh,vi,xyr]= simulate_2d_movement(32, smdl );
% Only keep 4 to 14 (rest are too far away)
idx= 4:14;
vi= vi(:,idx); xyr= xyr(:,idx);

clf;subplot(121)
show_fem(smdl); axis square

% Show target positions
   theta= linspace(0,2*pi,50); xr= cos(theta); yr= sin(theta);
   hold on;
   for i=1:size(xyr,2)
       hh= plot(xyr(3,i)*xr+ xyr(1,i),xyr(3,i)*yr- xyr(2,i));
       set(hh,'LineWidth',3,'Color',[0,0,0]);
   end
   hold off;



print -r125 -dpng dual_partial2d02a.png;
