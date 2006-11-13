% Lung images
% $Id: tutorial310b.m,v 1.1 2006-11-13 18:03:37 aadler Exp $

% If electrodes are counter-clockwise, then do this
imdl_ccw = imdl;
imdl_ccw.fwd_model.electrode([1,16:-1:2])=  ...
   imdl.fwd_model.electrode;

subplot(221);
show_fem(imdl_ccw.fwd_model);

% If electrodes start on back (dorsal), then do this
imdl_d = imdl;
imdl_d.fwd_model.electrode([9:16,1:8])=  ...
   imdl.fwd_model.electrode;

subplot(222);
show_fem(imdl_d.fwd_model);

axis equal
print -r100 -dpng tutorial310b.png;
