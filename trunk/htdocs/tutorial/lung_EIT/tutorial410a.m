% Lung images
% $Id: tutorial410a.m,v 1.1 2008-03-13 21:34:55 aadler Exp $

% 2D Model
imdl= mk_common_model('c2t3',16);


% most EIT systems image best with normalized difference
imdl.fwd_model.normalized_measurements= 1;
imdl.RtR_prior= @gaussian_HPF_prior;

% electrodes start on back (dorsal), then do this
imdl.fwd_model.electrode([9:16,1:8])=  ...
   imdl.fwd_model.electrode;


subplot(221);
show_fem(imdl.fwd_model);

axis equal
print -r100 -dpng tutorial410a.png;
